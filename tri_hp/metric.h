//
//  metric.h
//  tri_hp
//
//  Created by Brian Helenbrook on 3/2/25.
//

#ifndef metric_h
#define metric_h

#include <mapped_mesh.h>
#include "tri_hp.h"

class mapped_metric : public tri_hp::metric {
public:
    shared_ptr<mapping> map;
    mapped_metric(tri_hp& xin) : metric(xin) {}
    mapped_metric(const mapped_metric& in_metric, tri_hp& xin) : metric(xin), map(in_metric.map) {}
    virtual std::unique_ptr<metric> create(tri_hp& xin) override {return std::make_unique<mapped_metric>(*this,xin);}
    virtual void init(input_map& inmap) override;
    virtual void calc_metrics(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND,tri_mesh::ND>& dcrd, int tlvl=0) const override;
    virtual void calc_metrics1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& dcrd, int tlvl=0) const override;
    virtual void calc_positions(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl=0) const override;
    virtual void calc_positions1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl=0) const override;
    virtual void calc_positions0D(int vind, TinyVector<FLT,tri_mesh::ND>& pt, int tlvl=0) const override;
    virtual void setinfo() override;
};

class allcurved_metric : public tri_hp::metric {
public:
    allcurved_metric(tri_hp& xin) : metric(xin) {}
    allcurved_metric(const allcurved_metric& in_metric, tri_hp& xin) : metric(xin) {}
    virtual std::unique_ptr<metric> create(tri_hp& xin) override {return std::make_unique<allcurved_metric>(*this,xin);}
    virtual void init(input_map& inmap) override;
    virtual void calc_metrics(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND,tri_mesh::ND>& dcrd, int tlvl=0) const override;
    virtual void calc_metrics1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& dcrd, int tlvl=0) const override;
    virtual void calc_positions(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl=0) const override;
    virtual void calc_positions1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl=0) const override;
    virtual void calc_positions0D(int vind, TinyVector<FLT,tri_mesh::ND>& pt, int tlvl=0) const override;
    virtual void setinfo() override;
};

#endif /* metric_h */
