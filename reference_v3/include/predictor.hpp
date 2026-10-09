#pragma once
#include "engine.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Experimental, explicit prefix-cost surrogate inspired by family 121.
// No correctness or complexity guarantee from that paper is assumed here.
// All prediction errors are handled by the independent engine certificate.
namespace cba {
namespace detail {

inline int block_distance(const std::string& x, int a, int b,
                          const std::string& y, int c, int d,
                          std::uint64_t& work) {
    const int n = b-a, m = d-c;
    std::vector<int> prev(m+1), curr(m+1);
    for (int j=0; j<=m; ++j) prev[j]=j;
    for (int i=1; i<=n; ++i) {
        curr[0]=i;
        for (int j=1; j<=m; ++j)
            curr[j]=std::min({prev[j]+1, curr[j-1]+1,
                             prev[j-1]+int(x[a+i-1]!=y[c+j-1])});
        prev.swap(curr);
    }
    work += std::uint64_t(n)*m;
    return prev[m];
}

inline std::vector<std::pair<int,int>> anchors(const std::string& x,
                                               const std::string& y,
                                               int k, int spacing) {
    const int n=int(x.size()), m=int(y.size());
    std::vector<std::pair<int,int>> out{{0,0}};
    if (k>n || k>m || n==0 || m==0) { out.emplace_back(n,m); return out; }
    // Unique exact seeds on both sequences; a monotone LIS is only a proposal,
    // never a constraint on the certified optimum.
    using Occ = std::pair<int,int>;
    std::unordered_map<std::string,Occ> ax, ay;
    for(int i=0;i+k<=n;++i) {
        auto ins=ax.emplace(x.substr(i,k),Occ{i,0}); ++ins.first->second.second;
    }
    for(int j=0;j+k<=m;++j) {
        auto ins=ay.emplace(y.substr(j,k),Occ{j,0}); ++ins.first->second.second;
    }
    std::vector<std::pair<int,int>> hits;
    for(int i=1;i+k<n;++i) {
        auto key=x.substr(i,k); auto a=ax.find(key); auto b=ay.find(key);
        if(a->second.second==1 && b!=ay.end() && b->second.second==1) {
            int j=b->second.first;
            if(j>0 && j<m) hits.emplace_back(i,j);
        }
    }
    std::vector<int> tails, indices, predecessor(hits.size(),-1);
    for(int i=0;i<int(hits.size());++i) {
        int p=int(std::lower_bound(tails.begin(),tails.end(),hits[i].second)-tails.begin());
        if(p>0) predecessor[i]=indices[p-1];
        if(p==int(tails.size())) { tails.push_back(hits[i].second); indices.push_back(i); }
        else { tails[p]=hits[i].second; indices[p]=i; }
    }
    std::vector<std::pair<int,int>> chain;
    if(!indices.empty())
        for(int i=indices.back();i>=0;i=predecessor[i]) chain.push_back(hits[i]);
    std::reverse(chain.begin(),chain.end());
    for(auto a:chain)
        if(a.first-out.back().first>=spacing && a.second-out.back().second>=spacing
           && n-a.first>=k && m-a.second>=k) out.push_back(a);
    out.emplace_back(n,m);
    return out;
}

inline void linear_segment(std::vector<double>& centers,int x0,int x1,int y0,int y1) {
    if(x1==x0) { centers[x0]=double(y1); return; }
    for(int i=x0;i<=x1;++i)
        centers[i]=double(y0)+double(i-x0)*double(y1-y0)/double(x1-x0);
}

inline void cost_segment(Prediction& pred,const std::string& x,const std::string& y,
                         int x0,int x1,int y0,int y1,int block) {
    const int a=x1-x0,b=y1-y0,delta=b-a;
    if(a==0 || b==0 || delta==0) {
        linear_segment(pred.centers,x0,x1,y0,y1); return;
    }
    const int count=(a+block-1)/block;
    std::vector<int> cuts(count+1), lengths(count);
    std::vector<double> mass(count), allocation(count);
    // Bound pilot work before running it. Falls back to geometry on large
    // rectangles; the certificate remains responsible for correctness.
    constexpr std::uint64_t pilot_cap=4000000;
    std::uint64_t estimated=0;
    for(int t=0;t<=count;++t) cuts[t]=std::min(a,t*block);
    for(int t=0;t<count;++t) {
        int l=cuts[t],r=cuts[t+1];
        int yl=int(std::int64_t(l)*b/a),yr=int(std::int64_t(r)*b/a);
        estimated+=std::uint64_t(r-l)*(yr-yl);
    }
    if(pred.probe_cells+estimated>pilot_cap) {
        linear_segment(pred.centers,x0,x1,y0,y1); return;
    }
    double total=0;
    for(int t=0;t<count;++t) {
        int l=cuts[t],r=cuts[t+1]; lengths[t]=r-l;
        int yl=int(std::int64_t(l)*b/a),yr=int(std::int64_t(r)*b/a);
        int cost=block_distance(x,x0+l,x0+r,y,y0+yl,y0+yr,pred.probe_cells);
        mass[t]=double(cost)+1.0; // fixed smoothing, independent of benchmark
        total+=mass[t];
    }
    if(delta>0) {
        for(int t=0;t<count;++t) allocation[t]=delta*mass[t]/total;
    } else {
        // Capacity-limited weighted deletion allocation preserves monotonicity.
        double lo=0,hi=1;
        auto amount=[&](double lambda) {
            double s=0; for(int t=0;t<count;++t) s+=std::min(double(lengths[t]),lambda*mass[t]);
            return s;
        };
        while(amount(hi)<-delta) hi*=2;
        for(int it=0;it<64;++it) {
            double mid=(lo+hi)/2; if(amount(mid)<-delta)lo=mid;else hi=mid;
        }
        for(int t=0;t<count;++t) allocation[t]=-std::min(double(lengths[t]),hi*mass[t]);
    }
    double target=y0;
    for(int t=0;t<count;++t) {
        double next=target+lengths[t]+allocation[t];
        if(t+1==count) next=y1;
        next=std::clamp(next,target,double(y1));
        for(int u=0;u<=lengths[t];++u)
            pred.centers[x0+cuts[t]+u]=target+(next-target)*u/lengths[t];
        target=next;
    }
}
} // namespace detail

inline Prediction predict(const std::string& x,const std::string& y,
                          const std::string& mode,int block,int seed_k) {
    if(block<=0 || block>4096 || seed_k<=0 || seed_k>128)
        throw std::invalid_argument("block must be 1..4096 and seed-k 1..128");
    if(mode!="full" && mode!="linear" && mode!="cost" &&
       mode!="anchor-linear" && mode!="anchor-cost")
        throw std::invalid_argument("unknown predictor mode");
    Prediction p;
    int n=int(x.size()),m=int(y.size()); p.centers.resize(n+1,0.0);
    std::vector<std::pair<int,int>> a{{0,0},{n,m}};
    if(mode=="anchor-linear" || mode=="anchor-cost")
        a=detail::anchors(x,y,seed_k,std::max(32,block*4));
    bool cost=mode=="cost" || mode=="anchor-cost";
    for(std::size_t z=1;z<a.size();++z) {
        if(cost) detail::cost_segment(p,x,y,a[z-1].first,a[z].first,a[z-1].second,a[z].second,block);
        else detail::linear_segment(p.centers,a[z-1].first,a[z].first,a[z-1].second,a[z].second);
    }
    if(n>0) {p.centers[0]=0;p.centers[n]=m;}
    return p;
}
} // namespace cba
