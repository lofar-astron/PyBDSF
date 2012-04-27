## ----------------------------------------------------------------------
## $Id$
## ----------------------------------------------------------------------

enable_language(Fortran)

# Always use -fPIC to enable linking with shared libs.
add_library(port3 STATIC
  dnsg.f dn2g.f drnsg.f drn2g.f
  d1mach.f da7sst.f dc7vfn.f dd7tpr.f dd7upd.f df7hes.f dg7lit.f dg7qts.f
  ditsum.f divset.f dl7itv.f dl7ivm.f dl7mst.f dl7nvr.f dl7sqr.f dl7srt.f
  dl7svn.f dl7svx.f dl7tsq.f dl7tvm.f dl7vml.f dn2cvp.f dn2lrd.f dn2rdp.f
  do7prd.f dparck.f dq7apl.f dq7rad.f dq7rfh.f dr7mdc.f drldst.f ds7cpr.f
  ds7lup.f ds7lvm.f dv2axy.f dv2nrm.f dv7cpy.f dv7dfl.f dv7prm.f dv7scl.f
  dv7scp.f dv7swp.f i1mach.f i7mdcn.f stopx.f)
set_target_properties(port3 PROPERTIES COMPILE_FLAGS "-fPIC")
