// params.h — tiny key = value input-file reader.
//
// Accepts AMReX-style `key = value` lines such as
//     ns.cfl = 0.5
//     geometry.prob_hi = 1.0 1.0
//     # comment
// Values are kept as strings and converted on request.  Unknown keys are
// ignored, so input files can carry comments and unused keys.
#pragma once
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iostream>

class Params {
public:
    // Parse `filename`; later duplicates of a key overwrite earlier ones.  Throws if the file cannot be opened.
    void read(const std::string& filename) {
        std::ifstream in(filename);
        if (!in) throw std::runtime_error("cannot open input file " + filename);
        std::string line;
        while (std::getline(in, line)) {
            auto hash = line.find('#'); if (hash != std::string::npos) line = line.substr(0, hash);
            auto eq = line.find('='); if (eq == std::string::npos) continue;
            std::string key = trim(line.substr(0, eq)), val = trim(line.substr(eq + 1));
            if (!key.empty()) kv_[key] = val;
        }
    }
    // True if key k was present in the input file.
    bool has(const std::string& k) const { return kv_.count(k) > 0; }
    // Value of key k converted to T, or `def` if absent.  Throws on a malformed value.
    template <class T> T get(const std::string& k, T def) const {
        auto it = kv_.find(k); if (it == kv_.end()) return def;
        std::istringstream ss(it->second); T v; ss >> v;
        if (ss.fail()) throw std::runtime_error("bad value for " + k + ": " + it->second);
        return v;
    }
    // String overload (so that get("key", "default") does not deduce const char*).
    std::string get(const std::string& k, const char* def) const {
        auto it = kv_.find(k); return it == kv_.end() ? std::string(def) : it->second;
    }
    // Whitespace-separated list of values, e.g. "1.0 1.0", or `def` if absent.
    template <class T> std::vector<T> get_vec(const std::string& k, std::vector<T> def) const {
        auto it = kv_.find(k); if (it == kv_.end()) return def;
        std::istringstream ss(it->second); std::vector<T> v; T x; while (ss >> x) v.push_back(x);
        return v;
    }
    // Echo all key = value pairs (sorted by key) to os.
    void print(std::ostream& os) const { for (auto& p : kv_) os << "  " << p.first << " = " << p.second << "\n"; }
private:
    // Strip leading/trailing whitespace.
    static std::string trim(const std::string& s) {
        auto a = s.find_first_not_of(" \t\r\n"); if (a == std::string::npos) return "";
        auto b = s.find_last_not_of(" \t\r\n"); return s.substr(a, b - a + 1);
    }
    std::map<std::string, std::string> kv_;
};
