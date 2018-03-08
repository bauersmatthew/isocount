/*
 * Author: Matthew Bauer
 */

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <regex>
#include <sys/stat.h>
#include <algorithm>
#include <iterator>
#include <queue>
#include <deque>

#include <tclap/CmdLine.h>

/* PROTOTYPES */

/** Split a string by a delimiter.
 *
 * If the last field is empty, it is ignored. */
template<char delim='\t'>
std::vector<std::string> split(std::string const& str);

/** Convert A to B using a stringstream. */
template<class A, class B>
B ss_conv(A const& a);

/** Convert the given string using a stringstream. */
template<class B>
B ss_conv_s(std::string const& a);

/** A flexible, trace-capable exception class. */
class MyErr {
public:
    /** Compose a message by concatenating arguments with a stringstream. */
    template<class... Args>
    MyErr(Args const&... args);

    /** Combine two errors in a stacktrace-like manner. */
    friend MyErr operator+(MyErr const& a, MyErr const& b);

    /** Print. */
    friend std::ostream& operator<<(std::ostream& os, MyErr const& err);

protected:
    std::string msg; /** The full error message. */
};

/** A structure to parse and store configuration info. */
struct Config {
    static std::string const VERSION; /** The program version. */

    /** Parse arguments from the command line. */
    Config(int argc, char **argv);

    std::string anno_path; /** The path to the BED6 annotation. */
    std::string alig_path; /** The path to the SAM alignment. */

    std::string defcut_spec; /** The cutoff description. */
    std::vector<std::string> cgrp_specs; /** Cutoff-group specifications. */

    std::string dist_dir; /** Whether to print the distribution. */
};

/** A convenience class representing a percent value. */
struct Percent {
    Percent(); /** 0%. */
    Percent(double d); /** Tests for validity. */

    operator double() const;
    friend std::istream& operator>>(std::istream& inp, Percent& p);

    double val;
};

/** Holds information about a genomic region. */
struct GenomeRegion {
    GenomeRegion(); /** {"", 0, 1}. */
    GenomeRegion(std::string const& chr, uint32_t point);
    GenomeRegion(std::string const& chr, uint32_t s, uint32_t e);

    /** Test if a is strictly less than b (no overlap).
     *
     * If a and b are on different references, the result is the alphabetical
     * ordering of the reference names. */
    friend bool operator<(GenomeRegion const& a, GenomeRegion const& b);

    std::string ref; /** The chromosome/contig/reference. */
    uint32_t start, end; /** Half-open interval: [start, end). */

protected:
    /** Validate the interval (throw on failure). */
    void validate() const;
};

/** One annotation feature. */
class Feature : public GenomeRegion {
public:
    /** Load from a BED6 record. */
    explicit Feature(std::string const& bed6rec);

    std::string name; /** Should be unique. */
};

/** An entire annotation.
 *
 * Use as a sorted vector of Features. */
class Annotation : public std::vector<Feature> {
public:
    /** Load from BED6. */
    explicit Annotation(std::string const& path);
};

/** Stores relevant information from a SAM alignment. */
class SAMRecord : public GenomeRegion {
public:
    SAMRecord();
    /** Parse from SAM text. */
    explicit SAMRecord(std::string const& rec);

    std::string cigar_raw; /** As given in the file. */

    using CigarOp = std::pair<uint32_t, char>;
    std::vector<CigarOp> cigar; /** Parsed into (count, op) pairs. */

    friend class SAMLoader;

private:
    operator bool() const; /* Test if the record is valid. */
    bool good;
};

/** Load a SAM file line by line. */
class SAMLoader {
public:
    explicit SAMLoader(std::string const& path);
    ~SAMLoader();

    SAMRecord next();

    /** Check if the file has ended. */
    operator bool() const;

protected:
    std::string fname;
    std::ifstream fin;
    int lnnum;
};

/** A 2-tailed cutoff. */
struct Cutoff {
    Cutoff();  /** 0, 0 */
    explicit Cutoff(std::string const& spec); /** Parse a cutoff spec. */

    Percent inclusion;
    Percent exclusion;
};

/** A table linking feature names to cutoffs. */
class CutoffTable : public std::unordered_map<std::string, Cutoff> {
public:
    CutoffTable(Annotation const& anno, std::string const& defspec,
            std::vector<std::string> const& cgrps);
};

/** Represents a splicing isoform. */
class Isoform {
public:
    Isoform(Annotation const& anno, SAMRecord const& rec);

    /** Get the score attached to each feature for this alignment.
     *
     * Score formula:
     * (matches - mismatches - deletions - insertions) / feature_length */
    std::unordered_map<std::string, double> get_scores() const;

protected:
    std::unordered_map<std::string, uint32_t> lengths;
    std::unordered_map<std::string, uint32_t> matches;
    std::unordered_map<std::string, uint32_t> mistakes;
};

/* IMPLEMENTATIONS */
template<char delim='\t'>
std::vector<std::string> split(std::string const& str) {
    std::vector<std::string> ret;
    std::string growing;
    for (char ch : str) {
        if (ch == delim) {
            ret.push_back(growing);
            growing.clear();
        }
        else {
            growing += ch;
        }
    }
    if (!growing.empty())
        ret.push_back(growing);
    return ret;
}

template<class A, class B>
B ss_conv(A const& a) {
    std::stringstream ss;
    B b;
    ss << a;
    ss >> b;
    if (!ss) {
        throw MyErr("failed to convert: ", a);
    }
    return b;
}
template<class B>
B ss_conv_s(std::string const& s) {
    return ss_conv<std::string, B>(s);
}

template<class... Args>
MyErr::MyErr(Args const&... args) {
    std::stringstream ss;
    using unpack = int[];
    unpack{0, ((void)(ss << args), 0)...};
    msg = ss.str();
}
MyErr operator+(MyErr const& a, MyErr const& b) {
    return MyErr(a.msg, "\n\t", b.msg);
}
std::ostream& operator<<(std::ostream& os, MyErr const& err) {
    os << "E: " << err.msg;
}

std::string const Config::VERSION = "1.0.0";
Config::Config(int argc, char **argv) {
    using namespace TCLAP;
    try {
        CmdLine cmd(
                "Quickly quantify the splice isoforms present in the alignment.",
                ' ', Config::VERSION);
        UnlabeledValueArg<std::string> anno_arg(
                "annotation", "The BED6 annotation.", true, "", "anno.bed");
        UnlabeledValueArg<std::string> alig_arg(
                "alignment", "The SAM alignment.", true, "", "alig.sam");
        ValueArg<std::string> defcut_arg(
                "c", "cutoff", "The default cutoff. Format: inclusion,exclusion",
                false, "", "i,e");
        MultiArg<std::string> cutgroups_arg(
                "g", "cutoff-group",
                "A cutoff group. Format: 'regex inclusion,exclusion'",
                false, "pat i,e");
        ValueArg<std::string> distdir_arg(
                "d", "dist",
                ("Output the distribution of mistake ratios for each feature "
                 "the given directory."),
                false, "", "dists/");
        cmd.add(anno_arg);
        cmd.add(alig_arg);
        cmd.add(defcut_arg);
        cmd.add(cutgroups_arg);
        cmd.add(distdir_arg);
        cmd.parse(argc, argv);

        anno_path = anno_arg.getValue();
        alig_path = alig_arg.getValue();
        defcut_spec = defcut_arg.getValue();
        cgrp_specs = cutgroups_arg.getValue();
        dist_dir = distdir_arg.getValue();

        if (dist_dir.empty() == defcut_spec.empty()) {
            // both or neither is specified
            throw MyErr("exactly one of -d or -c must be specified!");
        }
    }
    catch (ArgException const& err) {
        throw MyErr(err.error(), " for arg ", err.argId());
    }
}

Percent::Percent() : Percent(0.0) {}
Percent::Percent(double d) {
    if (d < 0.0 || d > 100.0) {
        throw MyErr("invalid percent: ", d);
    }
    val = d;
}
Percent::operator double() const {
    return val;
}
std::istream& operator>>(std::istream& inp, Percent& p) {
    inp >> p.val;
    if (p.val < 0.0 || p.val > 100.0) {
        throw MyErr("invalid percent: ", p.val);
    }
    return inp;
}

GenomeRegion::GenomeRegion() : GenomeRegion("", 0) {}
GenomeRegion::GenomeRegion(std::string const& r, uint32_t p)
    : GenomeRegion(r, p, p+1) {}
GenomeRegion::GenomeRegion(std::string const& r, uint32_t s, uint32_t e) {
    ref = r;
    start = s;
    end = e;
    validate();
}
bool operator<(GenomeRegion const& a, GenomeRegion const& b) {
    if (a.ref != b.ref) {
        return (a.ref < b.ref);
    }
    else {
        return a.end <= b.start;
    }
}
void GenomeRegion::validate() const {
    if (start >= end) {
        throw MyErr("invalid interval: (", start, ", ", end, ")");
    }
}

Feature::Feature(std::string const& bed6rec) {
    std::vector<std::string> fields = split<'\t'>(bed6rec);
    if (fields.size() < 4) {
        throw MyErr("incomplete BED6 record: ", bed6rec);
    }
    try {
        name = fields[3];

        ref = fields[0];
        start = ss_conv_s<uint32_t>(fields[1]);
        end = ss_conv_s<uint32_t>(fields[2]);
        validate();
    }
    catch (MyErr const& err) {
        throw err + MyErr("in parsing of BED6 feature record: ", bed6rec);
    }
}

Annotation::Annotation(std::string const& path) {
    std::ifstream fin(path);
    if (!fin) {
        throw MyErr("failed to open annotation: ", path);
    }

    std::string line;
    int lnnum = 0;
    try {
        while(lnnum++, std::getline(fin, line), fin) {
            if (!line.empty()) {
                emplace_back(line);
            }
        }
    }
    catch (MyErr const& err) {
        throw err + MyErr("in line ", lnnum, " of ", path);
    }

    if (empty()) {
        throw MyErr("annotation empty: ", path);
    }

    // validate
    std::sort(begin(), end());
    std::unordered_set<std::string> names;
    Feature *last = &(*(begin()));
    names.insert(last->name);
    for (auto it = std::next(begin()); it != end(); it++) {
        if (names.count(it->name)) {
            throw MyErr("feature name used more than once: ", it->name);
        }
        names.insert(it->name);

        if (last->end != it->start) {
            throw MyErr("annotation features not contiguous: ", last->name,
                    " and ", it->name);
        }
        last = &(*it);
    }
}

SAMRecord::SAMRecord() {
    good = false;
}
SAMRecord::SAMRecord(std::string const& rec) {
    if (rec.empty() || rec[0] == '@') { // ==> header line (or blank)
        good = false;
        return;
    }
    std::vector<std::string> fields = split<'\t'>(rec);
    if (rec.size() < 11) {
        throw MyErr("incomplete SAM record");
    }
    try {
        int flag = ss_conv_s<int>(fields[1]);
        if (flag & 0x4 || flag & 0x900) { // ==> unmapped or not primary
            good = false;
            return;
        }

        ref = fields[2];
        start = ss_conv_s<uint32_t>(fields[3]) - 1; // -1 b/c SAM coord is 1-based
        //end = start + ss_conv_s<uint32_t>(fields[8]);
        //validate();
        // (end calculation moved down)
    }
    catch (MyErr const& err) {
        throw err + MyErr("in parsing of SAM record");
    }

    cigar_raw = fields[5];
    // check that the cigar string contains = and X instead of M
    if (cigar_raw.find_first_of('M') != std::string::npos) {
        throw MyErr("SAM record does not contain match/mismatch info",
                "; was GMAP run with the --sam-extended-cigar option?");
    }
    
    // parse cigar
    std::string growing_count;
    uint32_t len;
    for (char ch : cigar_raw) {
        if (ch >= '0' && ch <= '9') {
            growing_count += ch;
        }
        else {
            uint32_t count = ss_conv_s<uint32_t>(growing_count);
            growing_count.clear();

            if (ch == 'D' || ch == 'N' || ch == '=' || ch == 'X') {
                len += count;
            }
            cigar.emplace_back(count, ch);
        }
    }

    // finally set the end here (for some reason the tlen field is not
    // always accurate?)
    end = start + len;
    validate();

    good = true;
}
SAMRecord::operator bool() const {
    return good;
}

SAMLoader::SAMLoader(std::string const& path) {
    fin.open(path);
    if (!fin) {
        throw MyErr("failed to open alignment file: ", path);
    }
    lnnum = 0;
    fname = path;
}
SAMLoader::~SAMLoader() {
    fin.close(); // this would be done automatically anyway but w/e
}
SAMRecord SAMLoader::next() {
    std::string line;
    SAMRecord rec;
    try {
        while (lnnum++, std::getline(fin, line), fin) {
            rec = SAMRecord(line);
            if (rec) {
                return rec;
            }
        }
    }
    catch (MyErr const& err) {
        throw err + MyErr("in line ", lnnum, " of ", fname);
    }
    return SAMRecord(); // returning here means we got to the end of the file
}
SAMLoader::operator bool() const {
    return (bool)fin;
}

Cutoff::Cutoff() {
    inclusion = 0;
    exclusion = 0;
}
Cutoff::Cutoff(std::string const& spec) {
    std::vector<std::string> fields = split<','>(spec);
    if (fields.size() != 2) {
        throw MyErr("incomplete cutoff spec: ", spec);
    }
    inclusion = ss_conv_s<Percent>(fields[0]);
    exclusion = ss_conv_s<Percent>(fields[1]);
    if (inclusion < exclusion) {
        throw MyErr("inclusion cutoff < exclusion cutoff in spec: ", spec);
    }
}

CutoffTable::CutoffTable(Annotation const& anno, std::string const& defspec,
        std::vector<std::string> const& cgrps) {
    Cutoff def(defspec);
    std::vector<std::pair<std::regex, Cutoff>> cgrps_compiled;
    for (std::string const& spec : cgrps) {
        std::vector<std::string> fields = split<' '>(spec);
        if (fields.size() != 2) {
            throw MyErr("incomplete cutoff group spec: '", spec, "'");
        }
        try {
            cgrps_compiled.emplace_back(std::regex(fields[0]), Cutoff(fields[1]));
        }
        catch (std::regex_error const& err) {
            throw MyErr("invalid regex pattern: ", fields[0])
                + MyErr("in parsing of cutoff group spec: '", spec, "'");
        }
        catch (MyErr const& err) {
            throw err + MyErr("in parsing of cutoff group spec: '", spec, "'");
        }
    }

    for (Feature const& f : anno) {
        // save cutoff for this feature
        Cutoff c = def;
        for (auto const& p : cgrps_compiled) {
            if (std::regex_match(f.name, p.first)) {
                c = p.second;
                break;
            }
        }
        emplace(f.name, c);
    }
}

Isoform::Isoform(Annotation const& anno, SAMRecord const& rec) {
    for (Feature const& f : anno) {
        matches[f.name] = 0;
        mistakes[f.name] = 0;
        lengths[f.name] = f.end-f.start;
    }

    uint32_t refpos = rec.start; // the current position in the reference
    
    // find the first feature containing this record
    auto upper = std::upper_bound(anno.begin(), anno.end(),
            GenomeRegion(rec.ref, refpos));
    auto lower = std::lower_bound(anno.begin(), anno.end(),
            GenomeRegion(rec.ref, refpos));
    if (lower == anno.end()) {
        // we are past everything! ==> will never match any features
        return;
    }
    // f_idx is the index of the first matched feature.
    // here, it must be in the range [-1, anno.size()-1];
    // -1        ==> before everything,
    // otherwise ==> refpos is in anno[f_idx]
    int f_idx = std::distance(anno.begin(), upper) -1;

    std::queue<uint32_t> feat_bounds; // relevant feature end-points
    if (f_idx == -1) {
        feat_bounds.push(anno[0].start);
    }
    for (int i = std::max(0, f_idx); i < anno.size(); i++) {
        feat_bounds.push(anno[i].end);
    }

    // a deque is more useful to us here (we may need to push ops back in)
    std::deque<SAMRecord::CigarOp> cigar(rec.cigar.begin(), rec.cigar.end());

    while (!cigar.empty() && f_idx < int(anno.size())) {
        uint32_t end = feat_bounds.front();
        feat_bounds.pop();

        uint32_t *match_count = nullptr;
        uint32_t *mistake_count = nullptr;
        if (f_idx != -1) {
            std::string const& fname = anno[f_idx].name;
            match_count = &(matches[fname]);
            mistake_count = &(mistakes[fname]);
        }

        // consume ops until we get to the end of this feature
        while (refpos < end && !cigar.empty()) {
            SAMRecord::CigarOp op = cigar.front();
            cigar.pop_front();

            uint32_t count = op.first;
            char ch = op.second;
            uint32_t consumed = 0; // consumed from SEQ (not ref)

            // the relevant tally (matches or mismatches)
            uint32_t *tally = nullptr;

            if (ch == 'I') {
                // ins; not ref-consuming
                consumed = count;
                tally = mistake_count;
            }
            else if (ch == 'D' || ch == 'X' || ch == 'N') {
                // del or mis; ref-consuming
                consumed = std::min(end-refpos, count);
                tally = mistake_count;
                refpos += consumed;
            }
            else if (ch == '=') {
                // match; ref-consuming
                consumed = std::min(end-refpos, count);
                tally = match_count;
                refpos += consumed;
            }
            else {
                // we don't care about anything else; just skip and move on
                consumed = count;
            }

            // add to the relevant tally (if f_idx != -1)
            if (tally) {
                *tally += consumed;
            }

            // push part of op back into the queue if we didn't consume all of
            // it on this feature.
            if (consumed != count) {
                cigar.emplace_front(count-consumed, ch);
            }
        }

        f_idx++;
    }
}
std::unordered_map<std::string, double> Isoform::get_scores() const {
    std::unordered_map<std::string, double> ret;
    for (auto const& p : lengths) {
        std::string const& fname = p.first;

        double f_len = p.second;
        double n_good = int(matches.at(fname));
        double n_bad = int(mistakes.at(fname));
        ret[fname] = (n_good-n_bad)/f_len;
    }
    return ret;
}

void mainsub_dists_only(Annotation const& anno, SAMLoader& alig,
        std::string const& ddir) {
    // table for storing distribution info
    std::unordered_map<std::string, std::vector<double>> dist_table;
    for (Feature const& f : anno) {
        dist_table[f.name] = std::vector<double>();
    }
    SAMRecord rec;
    while (rec = alig.next(), alig) {
        std::unordered_map<std::string, double> scores =
            Isoform(anno, rec).get_scores();
        for (auto const& p : scores) {
            dist_table[p.first].push_back(p.second);
        }
    }

    mkdir(ddir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // (fails if dir already exists)
    
    for (auto& p : dist_table) {
        std::string fpath = ddir + "/" + p.first + ".tsv";
        std::ofstream fout(fpath);
        if (!fout) {
            throw MyErr("failed to open distribution output file: ", fpath);
        }
        std::sort(p.second.begin(), p.second.end());
        for (double d : p.second) {
            fout << d << std::endl;
        }
        fout.close();
    }
}

int main(int argc, char **argv) {
    try {
        Config cfg(argc, argv);

        if (cfg.dist_dir.empty()) {
            throw MyErr("sorry -- isoform serialization is not yet supported.");
        }
        Annotation anno(cfg.anno_path);
        SAMLoader alig(cfg.alig_path);
        mainsub_dists_only(anno, alig, cfg.dist_dir);
    }
    catch (MyErr const& err) {
        std::cerr << err << std::endl;
        return -1;
    }
    return 0;
}
