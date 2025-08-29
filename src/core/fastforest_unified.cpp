/**

MIT License

Copyright (c) 2020 Jonas Rembser

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

#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <stdexcept>
#include <stdlib.h> /* strtol */
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

std::vector<TreeEnsembleResponseType> fastforest::FastForest::softmax(const FeatureType* array,
                                                                      TreeEnsembleResponseType baseResponse) const {
    std::vector<TreeEnsembleResponseType> out(nClasses());
    softmax(array, out.data(), baseResponse);
    return out;
}

void fastforest::FastForest::softmax(const FeatureType* array,
                                     TreeEnsembleResponseType* out,
                                     TreeEnsembleResponseType baseResponse) const {
    int nClass = nClasses();
    if (nClass <= 2) {
        throw std::runtime_error(
            "Error in FastForest::softmax : binary classification models don't support softmax evaluation. Please set "
            "the number of classes in the FastForest-creating function if this is a multiclassification model.");
    }

    evaluate(array, out, nClass, baseResponse);
    fastforest::details::softmaxTransformInplace(out, nClass);
}

void fastforest::FastForest::evaluate(const FeatureType* array,
                                      TreeEnsembleResponseType* out,
                                      int nOut,
                                      TreeEnsembleResponseType baseResponse) const {
    for (int i = 0; i < nOut; ++i) {
        out[i] = baseResponse + baseResponses_[i];
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

TreeEnsembleResponseType fastforest::FastForest::evaluateBinary(const FeatureType* array,
                                                                TreeEnsembleResponseType baseResponse) const {
    TreeEnsembleResponseType out = baseResponse + baseResponses_[0];

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

namespace {

    namespace util {

        inline bool isInteger(const std::string& s) {
            if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+')))
                return false;

            char* p;
            strtol(s.c_str(), &p, 10);

            return (*p == 0);
        }

        template <class NumericType>
        struct NumericAfterSubstrOutput {
            explicit NumericAfterSubstrOutput() {
                value = 0;
                found = false;
                failed = true;
            }
            NumericType value;
            bool found;
            bool failed;
            std::string rest;
        };

        template <class NumericType>
        inline NumericAfterSubstrOutput<NumericType> numericAfterSubstr(std::string const& str,
                                                                        std::string const& substr) {
            std::string rest;
            NumericAfterSubstrOutput<NumericType> output;
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

}  // namespace

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

    while (std::getline(file, line)) {
        std::size_t foundBegin = line.find("[");
        std::size_t foundEnd = line.find("]");
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
                    std::stringstream ss2(splitstring[1]);
                    ss2 >> cutValue;
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
                util::NumericAfterSubstrOutput<int> output = util::numericAfterSubstr<int>(line, "yes=");
                if (!output.failed) {
                    yes = output.value;
                } else {
                    throw std::runtime_error(info + "problem while parsing the text dump");
                }
                output = util::numericAfterSubstr<int>(output.rest, "no=");
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

        } else {
            util::NumericAfterSubstrOutput<TreeResponseType> output =
                util::numericAfterSubstr<TreeResponseType>(line, "leaf=");
            if (output.found) {
                std::stringstream ss(line);
                int index;
                ss >> index;
                line = ss.str();

                ff.responses_.push_back(output.value);
                std::size_t nLeafIndices = leafIndices.size();
                leafIndices[index] = nLeafIndices + nPreviousLeaves;
            }
        }
    }
    terminateTree(ff, nPreviousNodes, nPreviousLeaves, nodeIndices, leafIndices, treesSkipped);

    if (nClasses > 2 && (ff.rootIndices_.size() + treesSkipped) % nClasses != 0) {
        std::stringstream ss;
        ss << "Error in FastForest construction : Forest has " << ff.rootIndices_.size()
           << " trees, which is not compatible with " << nClasses << "classes!";
        throw std::runtime_error(ss.str());
    }

    return ff;
}
