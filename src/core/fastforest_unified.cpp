/**

MIT License

Copyright (c) 2025 Jonas Rembser

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "core/fastforest.hpp"

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <sstream>
#include <stdexcept>
#include <stdlib.h> /* strtol */
#include <streambuf>
#include <vector>

namespace fastforest {
    namespace detail {

        typedef std::map<int, int> IndexMap;

        void correctIndices(std::vector<int>::iterator begin,
                            std::vector<int>::iterator end,
                            IndexMap const& nodeIndices,
                            IndexMap const& leafIndices);

    }  // namespace detail

}  // namespace fastforest

void fastforest::detail::correctIndices(std::vector<int>::iterator begin,
                                        std::vector<int>::iterator end,
                                        fastforest::detail::IndexMap const& nodeIndices,
                                        fastforest::detail::IndexMap const& leafIndices) {
    for (std::vector<int>::iterator it = begin; it != end; ++it) {
        IndexMap::const_iterator foundNode = nodeIndices.find(*it);
        if (foundNode != nodeIndices.end()) {
            *it = foundNode->second;
            continue;
        }
        IndexMap::const_iterator foundLeaf = leafIndices.find(*it);
        if (foundLeaf != leafIndices.end()) {
            *it = -foundLeaf->second;
            continue;
        } else {
            throw std::runtime_error("something is wrong in the node structure");
        }
    }
}

using namespace fastforest;

void fastforest::details::softmaxTransformInplace(TreeEnsembleResponseType* out, int nOut) {
    // Do softmax transformation inplace, mimicking exactly the Softmax function
    // in the src/common/math.h source file of xgboost.
    double norm = 0.;
    TreeEnsembleResponseType wmax = *out;
    int i = 1;
    for (; i < nOut; ++i) {
        wmax = std::max(out[i], wmax);
    }
    i = 0;
    for (; i < nOut; ++i) {
        TreeEnsembleResponseType& x = out[i];
        x = std::exp(x - wmax);
        norm += x;
    }
    i = 0;
    for (; i < nOut; ++i) {
        out[i] /= static_cast<float>(norm);
    }
}

std::vector<TreeEnsembleResponseType> fastforest::FastForest::softmax(const FeatureType* array) const {
    std::vector<TreeEnsembleResponseType> out(nClasses());
    softmax(array, out.data());
    return out;
}

void fastforest::FastForest::softmax(const FeatureType* array, TreeEnsembleResponseType* out) const {
    int nClass = nClasses();
    if (nClass <= 2) {
        throw std::runtime_error(
            "Error in FastForest::softmax : binary classification models don't support softmax evaluation. Please set "
            "the number of classes in the FastForest-creating function if this is a multiclassification model.");
    }

    evaluate(array, out, nClass);
    fastforest::details::softmaxTransformInplace(out, nClass);
}

void fastforest::FastForest::evaluate(const FeatureType* array, TreeEnsembleResponseType* out, int nOut) const {
    for (int i = 0; i < nOut; ++i) {
        out[i] = baseResponses_[i];
    }

    int iRootIndex = 0;
    for (std::vector<int>::const_iterator indexIter = rootIndices_.begin(); indexIter != rootIndices_.end();
         ++indexIter) {
        int index = *indexIter;
        do {
            int r = rightIndices_[index];
            int l = leftIndices_[index];
            index = array[cutIndices_[index]] < cutValues_[index] ? l : r;
        } while (index > 0);
        out[treeNumbers_[iRootIndex] % nOut] += responses_[-index];
        ++iRootIndex;
    }
}

TreeEnsembleResponseType fastforest::FastForest::evaluateBinary(const FeatureType* array) const {
    TreeEnsembleResponseType out = baseResponses_[0];

    for (std::vector<int>::const_iterator indexIter = rootIndices_.begin(); indexIter != rootIndices_.end();
         ++indexIter) {
        int index = *indexIter;
        do {
            int r = rightIndices_[index];
            int l = leftIndices_[index];
            index = array[cutIndices_[index]] < cutValues_[index] ? l : r;
        } while (index > 0);
        out += responses_[-index];
    }

    return out;
}

FastForest fastforest::load_bin(std::string const& txtpath) {
    std::ifstream ifs(txtpath.c_str(), std::ios::binary);
    return load_bin(ifs);
}

FastForest fastforest::load_bin(std::istream& is) {
    FastForest ff;

    int nRootNodes;
    int nNodes;
    int nLeaves;

    is.read(reinterpret_cast<char*>(&nRootNodes), sizeof(int));
    is.read(reinterpret_cast<char*>(&nNodes), sizeof(int));
    is.read(reinterpret_cast<char*>(&nLeaves), sizeof(int));

    ff.rootIndices_.resize(nRootNodes);
    ff.cutIndices_.resize(nNodes);
    ff.cutValues_.resize(nNodes);
    ff.leftIndices_.resize(nNodes);
    ff.rightIndices_.resize(nNodes);
    ff.responses_.resize(nLeaves);
    ff.treeNumbers_.resize(nRootNodes);

    is.read(reinterpret_cast<char*>(ff.rootIndices_.data()), nRootNodes * sizeof(int));
    is.read(reinterpret_cast<char*>(ff.cutIndices_.data()), nNodes * sizeof(CutIndexType));
    is.read(reinterpret_cast<char*>(ff.cutValues_.data()), nNodes * sizeof(FeatureType));
    is.read(reinterpret_cast<char*>(ff.leftIndices_.data()), nNodes * sizeof(int));
    is.read(reinterpret_cast<char*>(ff.rightIndices_.data()), nNodes * sizeof(int));
    is.read(reinterpret_cast<char*>(ff.responses_.data()), nLeaves * sizeof(TreeResponseType));
    is.read(reinterpret_cast<char*>(ff.treeNumbers_.data()), nRootNodes * sizeof(int));

    int nBaseResponses;
    is.read(reinterpret_cast<char*>(&nBaseResponses), sizeof(int));
    ff.baseResponses_.resize(nBaseResponses);
    is.read(reinterpret_cast<char*>(ff.baseResponses_.data()), nBaseResponses * sizeof(TreeEnsembleResponseType));

    return ff;
}

void fastforest::FastForest::write_bin(std::string const& filename) const {
    std::ofstream os(filename.c_str(), std::ios::binary);

    int nRootNodes = rootIndices_.size();
    int nNodes = cutValues_.size();
    int nLeaves = responses_.size();
    int nBaseResponses = baseResponses_.size();

    os.write(reinterpret_cast<const char*>(&nRootNodes), sizeof(int));
    os.write(reinterpret_cast<const char*>(&nNodes), sizeof(int));
    os.write(reinterpret_cast<const char*>(&nLeaves), sizeof(int));
    os.write(reinterpret_cast<const char*>(rootIndices_.data()), nRootNodes * sizeof(int));
    os.write(reinterpret_cast<const char*>(cutIndices_.data()), nNodes * sizeof(CutIndexType));
    os.write(reinterpret_cast<const char*>(cutValues_.data()), nNodes * sizeof(FeatureType));
    os.write(reinterpret_cast<const char*>(leftIndices_.data()), nNodes * sizeof(int));
    os.write(reinterpret_cast<const char*>(rightIndices_.data()), nNodes * sizeof(int));
    os.write(reinterpret_cast<const char*>(responses_.data()), nLeaves * sizeof(TreeResponseType));
    os.write(reinterpret_cast<const char*>(treeNumbers_.data()), nRootNodes * sizeof(int));
    os.write(reinterpret_cast<const char*>(&nBaseResponses), sizeof(int));
    os.write(reinterpret_cast<const char*>(baseResponses_.data()), nBaseResponses * sizeof(TreeEnsembleResponseType));
    os.close();
}

namespace util {

    // ------------------- float version -------------------
    float nextafter(float x, float y) {
        if (x != x || y != y)
            return std::numeric_limits<float>::quiet_NaN();
        if (x == y)
            return y;
        if (x == 0.0f) {
            float smallest = std::numeric_limits<float>::denorm_min();
            return (y > 0) ? smallest : -smallest;
        }

        int32_t bits;
        std::memcpy(&bits, &x, sizeof(float));

        if ((x > 0) == (y > x)) {
            ++bits;
        } else {
            --bits;
        }

        float result;
        std::memcpy(&result, &bits, sizeof(float));
        return result;
    }

    // ------------------- double version -------------------
    //double nextafter(double x, double y) {
    //    if (x != x || y != y)
    //        return std::numeric_limits<double>::quiet_NaN();
    //    if (x == y)
    //        return y;
    //    if (x == 0.0) {
    //        double smallest = std::numeric_limits<double>::denorm_min();
    //        return (y > 0) ? smallest : -smallest;
    //    }

    //    int64_t bits;
    //    std::memcpy(&bits, &x, sizeof(double));

    //    if ((x > 0) == (y > x)) {
    //        ++bits;
    //    } else {
    //        --bits;
    //    }

    //    double result;
    //    std::memcpy(&result, &bits, sizeof(double));
    //    return result;
    //}

    inline bool isInteger(const std::string& s) {
        if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+')))
            return false;

        char* p;
        strtol(s.c_str(), &p, 10);

        return (*p == 0);
    }

    template <class Type_t>
    struct AfterSubstrOutput {
        explicit AfterSubstrOutput() {
            value = Type_t();
            found = false;
            failed = true;
        }
        Type_t value;
        bool found;
        bool failed;
        std::string rest;
    };

    template <class Type_t>
    inline AfterSubstrOutput<Type_t> afterSubstr(std::string const& str, std::string const& substr) {
        std::string rest;
        AfterSubstrOutput<Type_t> output;
        output.rest = str;

        std::size_t found = str.find(substr);
        if (found != std::string::npos) {
            output.found = true;
            std::stringstream ss(str.substr(found + substr.size(), str.size() - found + substr.size()));
            ss >> output.value;
            if (!ss.fail()) {
                output.failed = false;
                output.rest = ss.str();
            }
        }
        return output;
    }

    std::vector<std::string> split(std::string const& strToSplit, char delimiter) {
        std::stringstream ss(strToSplit);
        std::string item;
        std::vector<std::string> splitStrings;
        while (std::getline(ss, item, delimiter)) {
            splitStrings.push_back(item);
        }
        return splitStrings;
    }

    bool exists(std::string const& filename) {
        if (FILE* file = fopen(filename.c_str(), "r")) {
            fclose(file);
            return true;
        } else {
            return false;
        }
    }

}  // namespace util

void terminateTree(fastforest::FastForest& ff,
                   int& nPreviousNodes,
                   int& nPreviousLeaves,
                   fastforest::detail::IndexMap& nodeIndices,
                   fastforest::detail::IndexMap& leafIndices,
                   int& treesSkipped) {
    using namespace fastforest::detail;
    correctIndices(ff.rightIndices_.begin() + nPreviousNodes, ff.rightIndices_.end(), nodeIndices, leafIndices);
    correctIndices(ff.leftIndices_.begin() + nPreviousNodes, ff.leftIndices_.end(), nodeIndices, leafIndices);

    if (nPreviousNodes != static_cast<int>(ff.cutValues_.size())) {
        ff.treeNumbers_.push_back(ff.rootIndices_.size() + treesSkipped);
        ff.rootIndices_.push_back(nPreviousNodes);
    } else {
        int treeNumbers = ff.rootIndices_.size() + treesSkipped;
        ++treesSkipped;
        ff.baseResponses_[treeNumbers % ff.baseResponses_.size()] += ff.responses_.back();
        ff.responses_.pop_back();
    }

    nodeIndices.clear();
    leafIndices.clear();
    nPreviousNodes = ff.cutValues_.size();
    nPreviousLeaves = ff.responses_.size();
}

FastForest fastforest::load_txt(std::string const& txtpath, std::vector<std::string>& features, int nClasses) {
    const std::string info = "constructing FastForest from " + txtpath + ": ";

    if (!util::exists(txtpath)) {
        throw std::runtime_error(info + "file does not exists");
    }

    std::ifstream file(txtpath.c_str());
    return load_txt(file, features, nClasses);
}

FastForest fastforest::load_txt(std::istream& file, std::vector<std::string>& features, int nClasses) {
    if (nClasses < 2) {
        throw std::runtime_error("Error in fastforest::load_txt : nClasses has to be at least two");
    }

    const std::string info = "constructing FastForest from istream: ";

    FastForest ff;
    ff.baseResponses_.resize(nClasses == 2 ? 1 : nClasses);

    int treesSkipped = 0;

    int nVariables = 0;
    std::map<std::string, int> varIndices;
    bool fixFeatures = false;

    if (!features.empty()) {
        fixFeatures = true;
        nVariables = features.size();
        for (int i = 0; i < nVariables; ++i) {
            varIndices[features[i]] = i;
        }
    }

    std::string line;

    fastforest::detail::IndexMap nodeIndices;
    fastforest::detail::IndexMap leafIndices;

    int nPreviousNodes = 0;
    int nPreviousLeaves = 0;

    bool baseScoreDefined = false;
    TreeEnsembleResponseType baseScore = 0.0;

    while (std::getline(file, line)) {
        std::size_t foundBegin = line.find("[");
        std::size_t foundEnd = line.find("]");
        util::AfterSubstrOutput<TreeResponseType> leafOutput = util::afterSubstr<TreeResponseType>(line, "leaf=");
        util::AfterSubstrOutput<TreeEnsembleResponseType> baseScoreOutput =
            util::afterSubstr<TreeEnsembleResponseType>(line, "base_score=");
        if (foundBegin != std::string::npos) {
            std::string subline = line.substr(foundBegin + 1, foundEnd - foundBegin - 1);
            if (util::isInteger(subline) && !ff.responses_.empty()) {
                terminateTree(ff, nPreviousNodes, nPreviousLeaves, nodeIndices, leafIndices, treesSkipped);
            } else if (!util::isInteger(subline)) {
                std::stringstream ss(line);
                int index;
                ss >> index;
                line = ss.str();

                std::vector<std::string> splitstring = util::split(subline, '<');
                std::string const& varName = splitstring[0];
                FeatureType cutValue;
                {
                    bool lessEqual = false;
                    if (splitstring[1][0] == '=') {
                        splitstring[1] = splitstring[1].substr(1);
                        lessEqual = true;
                    }
                    ss = std::stringstream(splitstring[1]);
                    ss >> cutValue;
                    if (lessEqual) {
                        FeatureType val = cutValue;
                        cutValue = util::nextafter(val, std::numeric_limits<FeatureType>::infinity());
#if __cplusplus >= 201103L
                        assert(cutValue == std::nextafter(val, std::numeric_limits<FeatureType>::infinity()));
#endif
                    }
                }
                if (!varIndices.count(varName)) {
                    if (fixFeatures) {
                        throw std::runtime_error(info + "feature " + varName + " not in list of features");
                    }
                    varIndices[varName] = nVariables;
                    features.push_back(varName);
                    ++nVariables;
                }
                int yes;
                int no;
                util::AfterSubstrOutput<int> output = util::afterSubstr<int>(line, "yes=");
                if (!output.failed) {
                    yes = output.value;
                } else {
                    throw std::runtime_error(info + "problem while parsing the text dump");
                }
                output = util::afterSubstr<int>(output.rest, "no=");
                if (!output.failed) {
                    no = output.value;
                } else {
                    throw std::runtime_error(info + "problem while parsing the text dump");
                }

                ff.cutValues_.push_back(cutValue);
                ff.cutIndices_.push_back(varIndices[varName]);
                ff.leftIndices_.push_back(yes);
                ff.rightIndices_.push_back(no);
                std::size_t nNodeIndices = nodeIndices.size();
                nodeIndices[index] = nNodeIndices + nPreviousNodes;
            }

        } else if (leafOutput.found) {
            std::stringstream ss(line);
            int index;
            ss >> index;
            line = ss.str();

            ff.responses_.push_back(leafOutput.value);
            std::size_t nLeafIndices = leafIndices.size();
            leafIndices[index] = nLeafIndices + nPreviousLeaves;
        } else if (baseScoreOutput.found) {
            baseScore = baseScoreOutput.value;
            baseScoreDefined = true;
        }
    }
    terminateTree(ff, nPreviousNodes, nPreviousLeaves, nodeIndices, leafIndices, treesSkipped);

    if (!baseScoreDefined) {
        std::stringstream ss;
        ss << "\nERROR: The model dump is missing the required base_score=<float> line.\n"
           << "       Without this hint, FastForest cannot guarantee correct parsing,\n"
           << "       and inference results may be silently incorrect.\n\n"
           << "To ensure the version hint is always consistent with the actual model\n"
           << "you trained, we recommend appending the required line\n"
           << "right after dumping the model. For example:\n\n"
           << "    outfile = \"model.txt\"\n"
           << "    booster = model.get_booster()\n"
           << "    # Dump the model to a .txt file\n"
           << "    booster.dump_model(outfile, fmap=\"\", with_stats=False, dump_format=\"text\")\n"
           << "    # Append the base score (unfortunately missing in the .txt dump)\n"
           << "    with open(outfile, \"a\") as f:\n"
           << "        import json\n"
           << "\n"
           << "        json_dump = json.loads(booster.save_config())\n"
           << "        base_score = json_dump[\"learner\"][\"learner_model_param\"][\"base_score\"]\n"
           << "        f.write(f\"base_score={base_score}\\n\")\n\n";
        throw std::runtime_error(ss.str());
    }

    if (nClasses > 2 && (ff.rootIndices_.size() + treesSkipped) % nClasses != 0) {
        std::stringstream ss;
        ss << "Error in FastForest construction : Forest has " << ff.rootIndices_.size()
           << " trees, which is not compatible with " << nClasses << "classes!";
        throw std::runtime_error(ss.str());
    }

    for (std::size_t i = 0; i < ff.baseResponses_.size(); ++i) {
        ff.baseResponses_[i] += baseScore;
    }

    return ff;
}

namespace util {

    std::string readFile(const char* filename);

}

std::string util::readFile(const char* filename) {
    std::ifstream t(filename);
    std::string str;

    t.seekg(0, std::ios::end);
    str.reserve(t.tellg());
    t.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

    return str;
}

namespace tmva {
    class XMLAttributes {
      public:
        XMLAttributes() {
            boostWeight_ = NULL;
            itree_ = NULL;
            pos_ = NULL;
            depth_ = NULL;
            IVar_ = NULL;
            Cut_ = NULL;
            res_ = NULL;
            nType_ = NULL;
        };

        XMLAttributes(XMLAttributes const& other) {
            boostWeight_ = NULL;
            itree_ = NULL;
            pos_ = NULL;
            depth_ = NULL;
            IVar_ = NULL;
            Cut_ = NULL;
            res_ = NULL;
            nType_ = NULL;
            *this = other;
        };

        XMLAttributes& operator=(XMLAttributes const& other) {
            if (other.boostWeight_ != NULL) {
                boostWeight_ = new double;
                *boostWeight_ = *other.boostWeight_;
            }
            if (other.itree_ != NULL) {
                itree_ = new int;
                *itree_ = *other.itree_;
            }
            if (other.pos_ != NULL) {
                pos_ = new char;
                *pos_ = *other.pos_;
            }
            if (other.depth_ != NULL) {
                depth_ = new int;
                *depth_ = *other.depth_;
            }
            if (other.IVar_ != NULL) {
                IVar_ = new int;
                *IVar_ = *other.IVar_;
            }
            if (other.Cut_ != NULL) {
                Cut_ = new double;
                *Cut_ = *other.Cut_;
            }
            if (other.res_ != NULL) {
                res_ = new double;
                *res_ = *other.res_;
            }
            if (other.nType_ != NULL) {
                nType_ = new int;
                *nType_ = *other.nType_;
            }
            return *this;
        }

        ~XMLAttributes() { reset(); }

        // If we set an attribute that is already set, this will do nothing and return false.
        // Therefore an attribute has repeated and we know a new node has started.
        void set(std::string const& name, std::string const& value) {
            if (name == "itree")
                return setValue(itree_, value);
            if (name == "boostWeight")
                return setValue(boostWeight_, value);
            if (name == "pos")
                return setValue(pos_, value.substr(0, 1));
            if (name == "depth")
                return setValue(depth_, value);
            if (name == "IVar")
                return setValue(IVar_, value);
            if (name == "Cut")
                return setValue(Cut_, value);
            if (name == "res")
                return setValue(res_, value);
            if (name == "nType")
                return setValue(nType_, value);
        }

        bool hasValue(std::string const& name) {
            if (name == "itree")
                return itree_ != NULL;
            if (name == "boostWeight")
                return boostWeight_ != NULL;
            if (name == "pos")
                return pos_ != NULL;
            if (name == "depth")
                return depth_ != NULL;
            if (name == "IVar")
                return IVar_ != NULL;
            if (name == "Cut")
                return Cut_ != NULL;
            if (name == "res")
                return res_ != NULL;
            if (name == "nType")
                return nType_ != NULL;
            return false;
        }

        int const* itree() const { return itree_; };
        double const* boostWeight() const { return boostWeight_; };
        char const* pos() const { return pos_; };
        int const* depth() const { return depth_; };
        int const* IVar() const { return IVar_; };
        double const* Cut() const { return Cut_; };
        double const* res() const { return res_; };
        int const* nType() const { return nType_; };

        void reset() {
            if (boostWeight_ != NULL) {
                delete boostWeight_;
                boostWeight_ = NULL;
            }
            if (itree_ != NULL) {
                delete itree_;
                itree_ = NULL;
            }
            if (pos_ != NULL) {
                delete pos_;
                pos_ = NULL;
            }
            if (depth_ != NULL) {
                delete depth_;
                depth_ = NULL;
            }
            if (IVar_ != NULL) {
                delete IVar_;
                IVar_ = NULL;
            }
            if (Cut_ != NULL) {
                delete Cut_;
                Cut_ = NULL;
            }
            if (res_ != NULL) {
                delete res_;
                res_ = NULL;
            }
            if (nType_ != NULL) {
                delete nType_;
                nType_ = NULL;
            }
        }

      private:
        template <class T>
        void setValue(T*& member, std::string const& value) {
            if (!member) {
                member = new T;
            }
            std::stringstream ss(value);
            ss >> *member;
        }

        // from the tree root node node
        double* boostWeight_;
        int* itree_;
        char* pos_;
        int* depth_;
        int* IVar_;
        double* Cut_;
        double* res_;
        int* nType_;
    };

    struct BDTWithXMLAttributes {
        std::vector<double> boostWeights;
        std::vector<std::vector<XMLAttributes> > nodes;
    };

    BDTWithXMLAttributes readXMLFile(std::string const& filename);

    BDTWithXMLAttributes readXMLFile(std::string const& filename) {
        const std::string str = util::readFile(filename.c_str());

        std::size_t pos1 = 0;

        std::string name;
        std::string value;

        BDTWithXMLAttributes bdtXmlAttributes;

        std::vector<XMLAttributes>* currentTree = NULL;

        XMLAttributes* attrs = NULL;

        while ((pos1 = str.find('=', pos1)) != std::string::npos) {
            std::size_t pos2 = str.rfind(' ', pos1) + 1;

            name = str.substr(pos2, pos1 - pos2);

            pos2 = pos1 + 2;
            pos1 = str.find('"', pos2);

            value = str.substr(pos2, pos1 - pos2);

            if (name == "boostWeight") {
                bdtXmlAttributes.boostWeights.push_back(0.0);
                std::stringstream(value) >> bdtXmlAttributes.boostWeights.back();
            }

            if (name == "itree") {
                bdtXmlAttributes.nodes.push_back(std::vector<tmva::XMLAttributes>());
                currentTree = &bdtXmlAttributes.nodes.back();
                currentTree->push_back(tmva::XMLAttributes());
                attrs = &currentTree->back();
            }

            if (bdtXmlAttributes.nodes.empty())
                continue;

            if (attrs->hasValue(name)) {
                currentTree->push_back(tmva::XMLAttributes());
                attrs = &currentTree->back();
            }

            attrs->set(name, value);
        }

        if (bdtXmlAttributes.nodes.size() != bdtXmlAttributes.boostWeights.size()) {
            throw std::runtime_error("nodes size and bosstWeights size don't match");
        }

        return bdtXmlAttributes;
    }

}  // namespace tmva

using namespace fastforest;

namespace {

    struct SlowTreeNode {
        SlowTreeNode() {
            isLeaf = false;
            depth = -1;
            index = -1;
            yes = -1;
            no = -1;
            missing = -1;
            cutIndex = -1;
            cutValue = 0.0;
            leafValue = 0.0;
        }
        bool isLeaf;
        int depth;
        int index;
        int yes;
        int no;
        int missing;
        int cutIndex;
        double cutValue;
        double leafValue;
    };

    std::vector<SlowTreeNode> getSlowTreeNodes(std::vector<tmva::XMLAttributes> const& nodes) {
        std::vector<SlowTreeNode> xgbNodes(nodes.size());

        int xgbIndex = 0;
        for (int depth = 0; xgbIndex != static_cast<int>(nodes.size()); ++depth) {
            int iNode = 0;
            for (std::vector<tmva::XMLAttributes>::const_iterator node = nodes.begin(); node != nodes.end(); ++node) {
                if (*node->depth() == depth) {
                    xgbNodes[iNode].index = xgbIndex;
                    ++xgbIndex;
                }
                ++iNode;
            }
        }

        int iNode = 0;
        for (std::vector<tmva::XMLAttributes>::const_iterator node = nodes.begin(); node != nodes.end(); ++node) {
            SlowTreeNode& xgbNode = xgbNodes[iNode];
            xgbNode.isLeaf = *node->nType() != 0;
            xgbNode.depth = *node->depth();
            xgbNode.cutIndex = *node->IVar();
            xgbNode.cutValue = *node->Cut();
            xgbNode.leafValue = *node->res();
            if (!xgbNode.isLeaf) {
                xgbNode.yes = xgbNodes[iNode + 1].index;
                xgbNode.no = xgbNode.yes + 1;
                xgbNode.missing = xgbNode.yes;
            }
            ++iNode;
        }

        return xgbNodes;
    }

    typedef std::vector<SlowTreeNode> SlowTree;
    typedef std::vector<SlowTree> SlowForest;

    //std::ostream& operator<<(std::ostream& os, SlowTreeNode const& node) {
    //    for (int i = 0; i < node.depth; ++i) {
    //        os << "\t";
    //    }
    //    if (node.isLeaf) {
    //        os << node.index << ":leaf=" << node.leafValue;
    //    } else {
    //        os << node.index << ":[f" << node.cutIndex << "<" << node.cutValue << "]";
    //        os << " yes=" << node.yes << ",no=" << node.no << ",missing=" << node.missing;
    //    }
    //    return os;
    //}

    //std::ostream& operator<<(std::ostream& os, SlowTree const& nodes) {
    //    for (SlowTree::const_iterator node = nodes.begin(); node != nodes.end(); ++node) {
    //        os << *node << "\n";
    //    }
    //    return os;
    //}

    //std::ostream& operator<<(std::ostream& os, SlowForest const& forest) {
    //    int iTree = 0;
    //    for (SlowForest::const_iterator tree = forest.begin(); tree != forest.end(); ++tree) {
    //        os << "booster[" << iTree << "]:"
    //           << "\n";
    //        os << *tree;
    //        ++iTree;
    //    }
    //    return os;
    //}

}  // namespace

namespace fastforest {

    //FastForest load_slowforest(SlowForest const& xgb, std::vector<std::string>& features) {
    FastForest load_slowforest(SlowForest const& xgb) {
        FastForest ff;
        ff.baseResponses_.resize(2);

        //int nVariables = 0;
        std::map<std::string, int> varIndices;
        //bool fixFeatures = false;

        std::map<int, int> nodeIndices;
        std::map<int, int> leafIndices;

        int nPreviousNodes = 0;
        int nPreviousLeaves = 0;

        for (SlowForest::const_iterator tree = xgb.begin(); tree != xgb.end(); ++tree) {
            detail::correctIndices(
                ff.rightIndices_.begin() + nPreviousNodes, ff.rightIndices_.end(), nodeIndices, leafIndices);
            detail::correctIndices(
                ff.leftIndices_.begin() + nPreviousNodes, ff.leftIndices_.end(), nodeIndices, leafIndices);
            nodeIndices.clear();
            leafIndices.clear();
            nPreviousNodes = ff.cutValues_.size();
            nPreviousLeaves = ff.responses_.size();
            ff.rootIndices_.push_back(nPreviousNodes);
            for (SlowTree::const_iterator node = tree->begin(); node != tree->end(); ++node) {
                if (node->isLeaf) {
                    ff.responses_.push_back(node->leafValue);
                    leafIndices[node->index] = leafIndices.size() + nPreviousLeaves;
                } else {
                    ff.cutValues_.push_back(node->cutValue);
                    ff.cutIndices_.push_back(node->cutIndex);
                    ff.leftIndices_.push_back(node->yes);
                    ff.rightIndices_.push_back(node->no);
                    nodeIndices[node->index] = nodeIndices.size() + nPreviousNodes;
                }
            }
        }

        detail::correctIndices(
            ff.rightIndices_.begin() + nPreviousNodes, ff.rightIndices_.end(), nodeIndices, leafIndices);
        detail::correctIndices(
            ff.leftIndices_.begin() + nPreviousNodes, ff.leftIndices_.end(), nodeIndices, leafIndices);

        return ff;
    }

    //FastForest load_tmva_xml(std::string const& xmlpath, std::vector<std::string>& features) {
    FastForest load_tmva_xml(std::string const& xmlpath) {
        tmva::BDTWithXMLAttributes tmvaXML = tmva::readXMLFile(xmlpath);
        std::vector<std::vector<SlowTreeNode> > xgboostForest;
        for (std::vector<std::vector<tmva::XMLAttributes> >::const_iterator tree = tmvaXML.nodes.begin();
             tree != tmvaXML.nodes.end();
             ++tree) {
            xgboostForest.push_back(getSlowTreeNodes(*tree));
        }
        //return fastforest::load_slowforest(xgboostForest, features);
        return fastforest::load_slowforest(xgboostForest);
    }
}  // namespace fastforest
