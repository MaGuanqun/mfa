# example illustrating encoding/decoding with fixed number of control points

import diy
import mfa
import math

class PyBlock(mfa.Block):
    # debug
    def test(self, cp, caps=False):
        if caps:
            print("HELLO")
        else:
            print("hello")

# default program arguments
fun             = "sinc"
error           = True
dom_dim         = 1
pt_dim          = 2
geom_degree     = 1
vars_degree     = 3
ndom_pts        = 20
geom_nctrl_pts  = geom_degree + 1
vars_nctrl_pts  = 10

# default dataset arguments
d_args                  = mfa.DomainArgs(dom_dim, pt_dim)
d_args.weighted         = 0
d_args.n                = 0.0
d_args.multiblock       = False
d_args.verbose          = 1
# NB, arrays bound to STL vectors must be assigned wholesale, not modified elementwise
d_args.f                = [1.0]
d_args.geom_p           = [geom_degree]
d_args.vars_p           = [[vars_degree]]
d_args.ndom_pts         = [ndom_pts]
d_args.geom_nctrl_pts   = [geom_nctrl_pts]
d_args.vars_nctrl_pts   = [[vars_nctrl_pts]]
d_args.min              = [-4.0 * math.pi]
d_args.max              = [4.0 * math.pi]
d_args.s                = [10.0, 1.0]

# debug
# print(d_args)
# print(d_args.dom_dim, d_args.pt_dim, d_args.min, d_args.max, d_args.s)

# MPI, DIY world and master
w = diy.mpi.MPIComm()           # world
m = diy.Master(w)               # master

def add_block(gid, core, bounds, domain, link):
    b = PyBlock()
    b.init(core, domain, dom_dim, pt_dim, float(0.0))
    m.add(gid, b, link)

nblocks = w.size
domain = diy.ContinuousBounds(d_args.min, d_args.max)
d = diy.ContinuousDecomposer(dom_dim, domain, nblocks)
a = diy.ContiguousAssigner(w.size, nblocks)
d.decompose(w.rank, a, add_block)

# initialize input data
m.foreach(lambda b, cp: b.generate_analytical_data(cp, fun, d_args))

# compute the MFA
m.foreach(lambda b, cp: b.fixed_encode_block(cp, d_args))

# debug: compute error field
if error:
    m.foreach(lambda b, cp: b.range_error(cp, 1, True, True))

# print results
m.foreach(lambda b, cp: b.print_block(cp, True))

# save the results
diy.write_blocks("approx.out", m)

