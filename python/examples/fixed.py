# example illustrating encoding/decoding with fixed number of control points

import diy
import mfa
import math

# default program arguments
fun             = "sinc"
error           = True
dom_dim         = 2
pt_dim          = 3
geom_degree     = 1
vars_degree     = 4
ndom_pts        = 100
geom_nctrl_pts  = geom_degree + 1
vars_nctrl_pts  = 11

# default dataset arguments
d_args                  = mfa.DomainArgs(dom_dim, pt_dim)
d_args.weighted         = 0
d_args.n                = 0.0
d_args.multiblock       = False
d_args.verbose          = 1
# NB, arrays bound to STL vectors must be assigned wholesale, not modified elementwise
d_args.f                = [1.0, 1.0]
d_args.geom_p           = [geom_degree, geom_degree]
d_args.vars_p           = [[vars_degree, vars_degree]]
d_args.ndom_pts         = [ndom_pts, ndom_pts]
d_args.geom_nctrl_pts   = [geom_nctrl_pts, geom_nctrl_pts]
d_args.vars_nctrl_pts   = [[vars_nctrl_pts, vars_nctrl_pts]]
d_args.min              = [-4.0 * math.pi, -4.0 * math.pi]
d_args.max              = [4.0 * math.pi, 4.0 * math.pi]
d_args.s                = [10.0, 1.0, 1.0]

# MPI, DIY world and master
w = diy.mpi.MPIComm()           # world
m = diy.Master(w)               # master

# TODO: using python definition of add_block because b.add does not work
# add_block could be removed if b.add would work
def add_block(gid, core, bounds, domain, link):
    b = mfa.Block()
#     print(m)
#     mfa.Block.add(gid, core, bounds, domain, link, m, dom_dim, pt_dim, float(0.0))
#     print(gid, core, bounds, domain, link, m, dom_dim, pt_dim)
    b.init(core, domain, dom_dim, pt_dim, float(0.0))
    m.add(gid, b, link)

nblocks = w.size

# decompose domain using double precision bounds
domain = diy.DoubleContinuousBounds(d_args.min, d_args.max)
d = diy.DoubleContinuousDecomposer(dom_dim, domain, nblocks)
a = diy.ContiguousAssigner(w.size, nblocks)
# TODO: using python definition of add_block because b.add does not work
d.decompose(w.rank, a, add_block)
# TODO: following does not work
# d.decompose(w.rank, a, lambda gid, core, bounds, domain_, link: mfa.Block.add(gid, core, bounds, domain_, link, m, dom_dim, pt_dim, float(0.0)))

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
print("\n\nSaving blocks\n")
diy.write_blocks("approx.out", m, save = mfa.save_block)

# debug: load the results and print them out
# print("\n\nLoading blocks back in and printing them out\n")
# m1 = diy.Master(w)
# a1 = diy.ContiguousAssigner(w.size, -1)
# diy.read_blocks("approx.out", a1, m1, load = mfa.load_block)
# m1.foreach(lambda b,cp: b.print_block(cp, False))
