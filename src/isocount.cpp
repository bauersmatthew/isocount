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

    std::string defcut_spec; /** Default cutoff spec. */
    std::vector<std::string> cgrp_specs; /** Cutoff group specs. */

    bool antisense; /** Serialize isoforms in antisense order. */
    bool suppress_uninteresting; /** No (undecidable) or (empty) isoforms. */

    std::string dists_path; /** Directory to put distribution info in. */
    std::string isotbl_path; /** File to put the isoform table in. */
    std::string categ_path; /** File to put categorized read info in. */

    bool verbose_dists; /** Add read name to the distribution files. */
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

    std::string name; /** Query name. */
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
     * 100 * (matches / (matches + mistakes)) */
    std::unordered_map<std::string, double> get_scores() const;

    /** Use the given cutoff table to serialize (name) this isoform. */
    std::string serialize(CutoffTable const& cuts, bool reverse) const;

protected:
    struct FeatureInfo {
        FeatureInfo();
        uint32_t mat, mis, ins, del;
    };
    std::unordered_map<std::string, FeatureInfo> feat_info;
    Annotation const& anno;
};

/** Helper class for generating classic isoform table output. */
class IsoTableGen {
private:
    std::string out;
    CutoffTable const *cuttbl;
    bool rev, supp;

    std::unordered_map<std::string, uint32_t> iso_table;

public:
    IsoTableGen(std::string const& out, CutoffTable const *cuts,
            bool rev, bool suppress);
    void add(Isoform const& iso);
    void save();
};

/** Helper class for generating score distribution output. */
class ScoreDistGen {
private:
    std::string outdir;
    bool verbose;
    std::unordered_map<
        std::string,
        std::vector<std::pair<double, std::string>>> dist_table;

public:
    ScoreDistGen(std::string const& out, Annotation const& anno,
            bool verbose);
    void add(Isoform const& iso, SAMRecord const& rec);
    void save();
};

/** Helper class for generating read-specific isoform output. */
class CategTableGen {
private:
    std::string out;
    CutoffTable const *cuttbl;
    bool rev;
    std::vector<std::pair<std::string, std::string>> read_isos;
    // (read id --> iso name)
    
public:
    CategTableGen(std::string const& out, CutoffTable const *cuts,
            bool rev);
    void add(Isoform const& iso, SAMRecord const& rec);
    void save();
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

std::string const Config::VERSION = "1.1.0";
Config::Config(int argc, char **argv) {
    using namespace TCLAP;
    try {
        CmdLine cmd(
                "Quickly quantify the splice isoforms present in the alignment.",
                ' ', Config::VERSION);
        /* Inputs. */
        UnlabeledValueArg<std::string> anno_arg(
                "annotation", "The BED6 annotation.", true, "", "anno.bed");
        UnlabeledValueArg<std::string> alig_arg(
                "alignment", "The SAM alignment.", true, "", "alig.sam");
        cmd.add(anno_arg);
        cmd.add(alig_arg);

        /* Serialization cutoffs. */
        ValueArg<std::string> defcut_arg(
                "c", "cutoff", "The default cutoff. Format: inclusion,exclusion",
                false, "", "i,e");
        MultiArg<std::string> cutgroups_arg(
                "g", "cutoff-group",
                "A cutoff group. Format: 'regex inclusion,exclusion'",
                false, "pat i,e");
        cmd.add(defcut_arg);
        cmd.add(cutgroups_arg);

        /* Isoform serialization settings. */
        SwitchArg antisense_arg(
                "r", "antisense",
                "Sort features in antisense order when serializing isoforms.");
        SwitchArg suppress_arg(
                "S", "suppress",
                "Do not report (undecidable) and (empty) serialized isoforms.");
        cmd.add(antisense_arg);
        cmd.add(suppress_arg);

        /* Distribution output settings. */
        SwitchArg verbosedists_arg(
                "D", "verbose-dists",
                ("Include the source read ID in a second column for each row "
                 "in the score distribution tables."));
        cmd.add(verbosedists_arg);

        ValueArg<std::string> distdir_arg(
                "d", "dist",
                ("Output the distribution of mistake ratios for each feature "
                 "the given directory."),
                false, "", "dists/");
        ValueArg<std::string> isotbl_arg(
                "o", "table-out",
                ("Output the isoform table to the given file. If not given, "
                 "but --cutoff is given, then /dev/stdout is used."),
                false, "", "isos.tsv");
        ValueArg<std::string> categ_arg(
                "C", "categorize",
                ("Categorize (serialize) each read, and output a table "
                 "linking read IDs to isoform serializations to the given "
                 "file."),
                false, "", "read-isos.tsv");
        cmd.add(distdir_arg);
        cmd.add(isotbl_arg);
        cmd.add(categ_arg);

        /* Parse and store. */
        cmd.parse(argc, argv);

        anno_path = anno_arg.getValue();
        alig_path = alig_arg.getValue();

        defcut_spec = defcut_arg.getValue();
        cgrp_specs = cutgroups_arg.getValue();

        antisense = antisense_arg.getValue();
        suppress_uninteresting = suppress_arg.getValue();

        verbose_dists = verbosedists_arg.getValue();

        dists_path = distdir_arg.getValue();
        isotbl_path = isotbl_arg.getValue();
        categ_path = categ_arg.getValue();

        /* Do some error checking. */
        if (isotbl_path.size() && defcut_spec.empty()) {
            throw MyErr("cannot serialize isoforms without a default cutoff!");
        }
        if (categ_path.size() && defcut_spec.empty()) {
            throw MyErr("cannot serialize isoforms without a default cutoff!");
        }
        if (defcut_spec.size() && isotbl_path.empty()) {
            isotbl_path = "/dev/stdout";
        }
        if (antisense && defcut_spec.empty()) {
            throw MyErr("-r requires serialization mode (-c)");
        }
        if (suppress_uninteresting && defcut_spec.empty()) {
            throw MyErr("-S requires serialization mode (-c)");
        }
        if (verbose_dists && dists_path.empty()) {
            throw MyErr("-D requires -d");
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
        name = fields[0];

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

Isoform::FeatureInfo::FeatureInfo() {
    mat = mis = ins = del = 0;
}
Isoform::Isoform(Annotation const& annot, SAMRecord const& rec)
    : anno(annot) {

    for (Feature const& f : anno) {
        feat_info.emplace(f.name, FeatureInfo());
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

    // handle each feature
    while (!cigar.empty() && f_idx < int(anno.size())) {
        uint32_t end = feat_bounds.front();
        feat_bounds.pop();

        // the currently-relevant info struct
        FeatureInfo *info = nullptr;
        if (f_idx != -1) {
            info = &(feat_info[anno[f_idx].name]);
        }

        // consume ops until we get to the end of this feature
        while (refpos < end && !cigar.empty()) {
            SAMRecord::CigarOp op = cigar.front();
            cigar.pop_front();

            uint32_t count = op.first;
            char ch = op.second;
            uint32_t consumed = 0; // consumed from SEQ (not ref)
            uint32_t delta = 0; // moved in ref

            // the relevant tally (matches or mismatches)
            uint32_t *tally = nullptr;

            consumed = std::min(end-refpos, count); // ONLY FOR REF-CONS OPS
            delta = consumed; // ONLY TRUE FOR REF-CONSUMING OPERATIONS!
            if (ch == '=') {
                // match; ref-consuming
                if (info) {
                    tally = &(info->mat);
                }
            }
            else if (ch == 'X') {
                // mismatch; ref-consuming
                if (info) {
                    tally = &(info->mis);
                }
            }
            else if (ch == 'I') {
                // insertion; NOT ref-consuming!!
                consumed = count;
                delta = 0;
                if (info) {
                    tally = &(info->ins);
                }
            }
            else if (ch == 'D' || ch == 'N') {
                // deletion; ref-consuming
                if (info) {
                    tally = &(info->del);
                }
            }
            else {
                // we don't care about anything else; just skip and move on
                consumed = count;
                delta = 0;
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

            refpos += delta;
        }

        if (cigar.empty() && info) {
            // the read ended, possibly in the middle of a feature
            // ==> count the uncovered part of the feature as a deletion
            info->del += end-refpos;
        }
        f_idx++;
    }
}
std::unordered_map<std::string, double> Isoform::get_scores() const {
    std::unordered_map<std::string, double> ret;
    for (auto const& p : feat_info) {
        std::string const& fname = p.first;
        FeatureInfo const& info = feat_info.at(fname);

        double n_good = info.mat;
        double n_bad = info.mis + info.ins + info.del;
        double denom = n_good+n_bad;

        if (denom < 0.9) { // denom == 0
            ret[fname] = 0;
        }
        else {
            ret[fname] = 100.0*(n_good/denom);
        }
    }
    return ret;
}
std::string Isoform::serialize(CutoffTable const& cuts, bool rev) const {
    auto scores = get_scores();
    std::vector<std::reference_wrapper<Feature const>> included;

    for (Feature const& f : anno) {
        double s = scores[f.name];
        Cutoff const& cut = cuts.at(f.name);
        if (s >= cut.inclusion) {
            included.emplace_back(f);
        }
        else if (s > cut.exclusion) {
            // it doesn't fit in either category; throw us out!
            return "(undecidable)";
        }
    }

    // sort by feature position
    auto sense_sorter = [](
            std::reference_wrapper<Feature const> a,
            std::reference_wrapper<Feature const> b) {
        return a < b;
    };
    auto antisense_sorter = [](
            std::reference_wrapper<Feature const> a,
            std::reference_wrapper<Feature const> b) {
        return b < a;
    };
    std::sort(included.begin(), included.end(),
            rev ? antisense_sorter : sense_sorter);

    // combine into string
    std::string ret;
    for (std::reference_wrapper<Feature const> const& f : included) {
        ret += f.get().name;
        ret += ",";
    }
    if (!ret.empty()) {
        // remove trailing comma
        ret = ret.substr(0, ret.size()-1);
    }
    else {
        ret = "(empty)";
    }
    return ret;
}

IsoTableGen::IsoTableGen(std::string const& o, CutoffTable const *c,
        bool r, bool s) : out(o), cuttbl(c), rev(r), supp(s) {}
void IsoTableGen::add(Isoform const& iso) {
    if (out.empty()) {
        return;
    }
    std::string ser = iso.serialize(*cuttbl, rev);
    if (!supp || (ser != "(undecidable)" && ser != "(empty)")) {
        if (!iso_table.count(ser)) {
            iso_table.emplace(ser, 0);
        }
        iso_table[ser]++;
    }
}
void IsoTableGen::save() {
    if (out.empty()) {
        return;
    }
    std::ofstream fout(out);
    if (!fout) {
        throw MyErr("failed to open output file: ", out);
    }
    std::vector<std::pair<std::string, uint32_t>> sorted(
            iso_table.begin(), iso_table.end());
    std::sort(sorted.begin(), sorted.end(),
            [](
                std::pair<std::string, double> const& a,
                std::pair<std::string, double> const& b)
            { return a.second > b.second; });
    for (auto const& p : sorted) {
        fout << p.first << '\t' << p.second << std::endl;
    }
    fout.close();
}

ScoreDistGen::ScoreDistGen(std::string const& o, Annotation const& a,
        bool v) : outdir(o), verbose(v) {
    if (outdir.empty()) {
        return;
    }
    for (Feature const& f : a) {
        dist_table.emplace(f.name,
                std::vector<std::pair<double, std::string>>());
    }
}
void ScoreDistGen::add(Isoform const& iso, SAMRecord const& rec) {
    if (outdir.empty()) {
        return;
    }
    auto scores = iso.get_scores();
    for (auto const& p : scores) {
        dist_table[p.first].emplace_back(p.second, rec.name);
    }
}
void ScoreDistGen::save() {
    if (outdir.empty()) {
        return;
    }
    mkdir(outdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // (fails if dir already exists)
    
    for (auto& p : dist_table) {
        std::string fpath = outdir + "/" + p.first + ".tsv";
        std::ofstream fout(fpath);
        if (!fout) {
            throw MyErr("failed to open distribution output file: ", fpath);
        }
        std::sort(p.second.begin(), p.second.end(),
                [](
                    std::pair<double, std::string> const& a,
                    std::pair<double, std::string> const& b)
                { return a.first < b.first; });
        for (auto const& p2 : p.second) {
            fout << p2.first;
            if (verbose) {
                fout << '\t' << p2.second;
            }
            fout << std::endl;
        }
        fout.close();
    }
}

CategTableGen::CategTableGen(std::string const& o, CutoffTable const *c,
        bool r) : out(o), cuttbl(c), rev(r) {}
void CategTableGen::add(Isoform const& iso, SAMRecord const& rec) {
    if (out.empty()) {
        return;
    }
    std::string ser = iso.serialize(*cuttbl, rev);
    read_isos.emplace_back(rec.name, ser);
}
void CategTableGen::save() {
    if (out.empty()) {
        return;
    }
    std::ofstream fout(out);
    if (!fout) {
        throw MyErr("failed to open output file: ", out);
    }
    for (auto const& p : read_isos) {
        fout << p.first << '\t' << p.second << std::endl;
    }
    fout.close();
}

int main(int argc, char **argv) {
    try {
        Config cfg(argc, argv);

        Annotation anno(cfg.anno_path);
        SAMLoader alig(cfg.alig_path);

        std::unique_ptr<CutoffTable> cutoffs(nullptr);
        if (!cfg.defcut_spec.empty()) {
            cutoffs.reset(new CutoffTable(anno, cfg.defcut_spec,
                        cfg.cgrp_specs));
        }

        IsoTableGen g_isotbl(cfg.isotbl_path, cutoffs.get(),
            cfg.antisense, cfg.suppress_uninteresting);
        ScoreDistGen g_dists(cfg.dists_path, anno, cfg.verbose_dists);
        CategTableGen g_categ(cfg.categ_path, cutoffs.get(), cfg.antisense);

        SAMRecord rec;
        while (rec = alig.next(), alig) {
            Isoform iso(anno, rec);

            g_isotbl.add(iso);
            g_dists.add(iso, rec);
            g_categ.add(iso, rec);
        }

        g_isotbl.save();
        g_dists.save();
        g_categ.save();
    }
    catch (MyErr const& err) {
        std::cerr << err << std::endl;
        return -1;
    }
    return 0;
}
