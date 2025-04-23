r"""
This code uses the method of S. Allen Broughton to compute the homology
of a group action on a Riemann surface given by a generating vector. 

This was version 0; the code does not at present have a way to make the 
matrices symplectic, and no optimisation has been carried out to improve
performance. At present the slow part (as determined via profiling) is finding
the kernel of M1.

This was version 1; still not symplectic, optimised the right kernel process 
partially. 

This is version 2; a large rewrite has occurred. We now work with cohomology
alone to get both the intersection form and the group action. 

This is version 3; the underlying simplices have been adjusted so that the code
can calculate the intersection form when t>3. 

This is version 4; the right inverse RI_2 no longer needs to be calculated as 
it has been determined theoretically. Moreover, the way of constructing the 
sparse RI_1 has been improved to reduce on overhead. 

AUTHORS:

- Linden Disney-Hogg (2024-07-24): initial working version

"""

def homology_representation(G, gv, BR=ZZ, intersection_matrix=False, elms=None,
                            check_consistency=False):
    r"""
    Return a homology representation determined by a generating vector, assumed
    to be acting with quotient genus zero, and genus >=2.

    INPUT:

    - ``G`` -- a permutation group.

    - ``gv`` -- a list containing elements of ``G`` which determine the 
      corresponding Riemann surface. These group elements should have trivial
      product and generate the group. 

    - ``BR`` -- (default: ``ZZ``). The ring over which to give the 
      representation. Must be able to map ``ZZ`` to this ring. 

    - ``intersection_matrix`` -- (default: False). Whether to compute the 
      intersection matrix of the cohomology basis used to compute the 
      representation. 

    - ``elms`` -- (default: ``None``). Which elements of the group to give the
      representation of. If ``None`` then the elements of ``gv`` are used. 

    - ``check_consistency`` -- (default: False). Whether to run checks as the 
      executes to search for bugs. This may slow down the method if set to 
      ``True`` when using large groups acting in high genus.

    OUTPUT:

    A list giving the homology representation of the elements of ``elms``. If
    ``intersection_matrix`` is set to ``True`` then an antisymmetric matrix
    giving the intersection form of the (co)homology basis is also returned.

    EXAMPLES:

        sage: G = CyclicPermutationGroup(2)
        sage: c = G.gen()
        sage: sigma = 2
        sage: gv = 2*(sigma + 1)*[c]
        sage: hr = homology_representation(G, gv, BR=ZZ, elms=[c])
        sage: hr[0] == -matrix.identity(2*sigma)
        True

    ALGORITHM: 

    We implement the algorithm of Broughton, see paper. 

    .. TODO::

        Investigate the use of GAP Simplicial Complex methodologies to speed up
        parts of the algorithm. More generally, using simiplicial complex data
        types that persistent homology researchers employ could speed up this 
        process. 

        At present the slowest part of the code for large groups is the
        computation of the right kernel matrix, but the PLUQ algorithm for
        dense matrices over GF(2) is VERY fast. If one can cast into this class
        efficiently then this could speed up the implementation. The problem
        is that naively asking for the dense version of the sparse matrix is 
        slow enough that any gains in calculating the right inverse are lost.
        There are also perhaps benefits to be found by using linbox, a C 
        package that historically was used in Sage to compute fast right
        kernels of sparse integer matrices, but was removed due to errors
        arising inside other parts of Sage as a result. 

        Use the fact that differentials are just built out of permutation 
        matrices to speed up their linear algebra. 
    
        A future improvement could consider how to make the code more generic
        to allow for arbitrary downstairs graphs. 
    """
    if elms is None:
        elms = gv

    # Initialise some properties of the signature. 
    GO = G.order()
    ns = [gj.order() for gj in gv]
    t = len(gv)
    sigma = 1 + ZZ(GO*(-2 + t - sum(1/nj for nj in ns))/2)
    # We will assume that sigma is >= 2 for this code. 
    assert sigma > 1

    # The following matrices are used to get the simplex bases from those that
    # are natural for the symmetry of the problem, to those that have the 
    # vertices correctly ordered within themselves. 
    IG = matrix.identity(BR, GO, sparse=True)
    # nde and ndf are shorthands for the number of downstairs edges and faces
    nde = t + 2*(t - 3)
    ndf = 2*(t - 2)
    edge_bc = block_diagonal_matrix([-IG] + (nde-1)*[IG])

    # To get M1, we realise that it is a H x |G|*nde matrix, acting as 
    # d1(gi*ej) = gi*vj - gi*v(j-1), where H is the sum of the coset sizes. 
    # The way that gi acts on vj is determined by the generating vector, and is
    # equivalent the action of left multiplication on the cosets.
    H = ZZ(GO*sum(1/ni for ni in ns))
    # https://docs.gap-system.org/doc/ref/chap41.html#X7FED50ED7ACA5FB2
    # These are the cosetactions, given a shorter name for simplicity. 
    # FactorCosetAction by default computes the action of right multiplication
    # on the right cosets Hg. To get the left action on the left cosets on must
    # take the inverse of permutations given by FactorCosetAction. The matrices
    # in Ps then store this projection from group elements to cosets. 
    G_gap = G.gap()
    fs = [G_gap.FactorCosetAction(G.subgroup([ci]).gap())
                    for ci in gv]
    Ps = [matrix.zero(BR, ZZ(GO/nj), GO, sparse=True) for nj in ns]
    for i, gi in enumerate(G):
        for fj, Pj in zip(fs, Ps):
            Pj[PermutationGroupElement(fj.ImageElm(gi)).inverse()(1)-1, i] = 1

    # We will also need the action of right multiplication by the hi to
    # construct the boundary matrices, so we create these as matrices rho which
    # act by left multiplication on vectors and store them.
    hs = [G.identity()]
    for j in range(t-1):
        hs.append(hs[-1]*gv[j])
    gs = G.list()
    rhos = [IG] + [matrix(BR, GO, GO, {(gs.index(gi*hj), i): 1 
                                for i, gi in enumerate(G)}, 
                          sparse=True)
                   for hj in hs[1:]]

    # We will want to build the differential M1 from blocks consisting of the 
    # projectors of group elements to cosets, including the hj multiplication
    # in the case that t>3 and we have had to add aditional edges from the cell
    # decomposition in order to triangulate. 
    # As the exact method for constructing these blocks will differ when t=3
    # and t>3, we insert a conditional and use the same condition for M2 also. 
    block_components = []
    for j in range(t):
        row_j = t*[0]
        row_j[j] = Ps[j]
        row_j[(j+1)%t] = -Ps[j]
        block_components.append(row_j)
    MM = block_matrix(BR, t, t, block_components, 
                      subdivide=False)
    if t > 3:
        # M1 part
        top_row = block_matrix(BR, 1, t-3, (t - 3)*[-Ps[0]],
                               subdivide=False)
        block_Ps = block_diagonal_matrix(Ps[2:-1])
        Prs = [Pj*rhoj for Pj, rhoj in zip(Ps, rhos)]
        block_Prs = block_diagonal_matrix(Prs[2:-1])
        Z2 = matrix.zero(BR, ZZ(GO/ns[1]), (t-3)*GO, sparse=True)
        Zt = matrix.zero(BR, ZZ(GO/ns[-1]), (t-3)*GO, sparse=True)
        M12 = block_matrix(BR, 4, 1, [[top_row], [Z2], [block_Ps], [Zt]], 
                           subdivide=False)
        M13 = block_matrix(BR, 4, 1, [[top_row], [Z2], [block_Prs], [Zt]], 
                           subdivide=False)
        M1 = block_matrix(BR, 1, 3, [MM, M12, M13])
        # M2 part
        Z = matrix.zero(BR, GO, (t-3)*GO, sparse=True)
        M21_top = block_matrix(2, 2, [[Z, IG], [IG, Z]],
                              subdivide=False)
        M22_top = block_matrix(2, 2, [[Z, rhos[0]], [rhos[1], Z]],
                              subdivide=False)
        block_rs = block_diagonal_matrix(rhos[2:])
        Id = matrix.identity(BR, (t-2)*GO, sparse=True)
        M21 = block_matrix(2, 1, [[M21_top], [Id]],
                          subdivide=False)
        M22 = block_matrix(2, 1, [[M22_top], [block_rs]],
                          subdivide=False)
        block_components = []
        for j in range(t-3):
            row_j = (t-2)*[0]
            row_j[j] = -IG
            row_j[(j+1)%(t-2)] = IG
            block_components.append(row_j)
        M23 = block_matrix(BR, t-3, t-2, block_components, 
                          subdivide=False)
        M2 = block_matrix(BR, 3, 2, [[M21, M22], [M23, 0], [0, M23]],
                          subdivide=False) 
    else:
        # M1 part 
        M1 = MM
        # M2 part
        M2 = block_matrix(3, 2, [[IG, rhos[0]],
                                 [IG, rhos[1]],
                                 [IG, rhos[2]]],
                         subdivide=False)
    # We transform to the new basis
    M1 = M1*edge_bc
    M2 = edge_bc*M2

    # We rename from M to N to make clear that we are now switching to a 
    # cohomology perspective. 
    N0, N1 = M1, M2
    
    # From Lemma 6 of SAB we know that the columns corresponding to the pivots
    # give a basis of the image. We construct the C part of the matrix R
    # using these 
    C = N0.matrix_from_rows(N0.pivot_rows())
    # Using right_kernel_matrix() instead of right_kernel().matrix() is 
    # considerably faster. 
    Dp = N1.transpose().right_kernel_matrix()
    Rp = block_matrix(2, 1, [C, Dp])
    R_cohom_1 = Rp.matrix_from_rows(Rp.pivot_rows())
    m_cohom_1 = C.nrows()
    D_cohom = R_cohom_1[m_cohom_1:, :]
    # We now calculate the appropriate partial right inverse. 
    # Efficiency can be gained here by shortcutting some Sage methods, and 
    # avoiding verifying the solution by assuming a known solution exists. 
    # It is likely future efficiency could be gained by not calculating the
    # full right inverse as an intermediate step. 
    nc = R_cohom_1.ncols()
    pivot_cols = R_cohom_1.pivots()
    # For future improvements, we don't need the full inverse here, just the 
    # latter columns. 
    X = R_cohom_1.matrix_from_columns(pivot_cols).inverse()
    if len(pivot_cols) < nc:
        # Because we are working with sparse matrices, this is the faster 
        # way to construct Y from the rows of X. 
        Xd = X.dict()
        Y = matrix(BR, nc, X.ncols(), {(pivot_cols[k[0]], k[1]): Xd[k]
                                       for k in Xd})
        RI_1 = Y[:, m_cohom_1:]

    # If we are going to compute the intersection matrix we need to pair up
    # elements of H^1 to make elements of H^2, and we need to read off the 
    # corresponding coefficient of this element in H^2 using a right inverse.  
    if intersection_matrix:
        # Define the right inverse of R_cohom_2 using known theorem.
        # This diagonal corresponds to the fundamental cycle. 
        RI_2_diag = matrix.diagonal(BR, (ndf*GO/2)*[1] + (ndf*GO/2)*[-1],
                                    sparse=True)
        # Calculate the IM
        # L_dict and R_dict give the behaviour of the map which restricts
        # 2-simplices to the left and right 1-simplices. 
        h2 = hs[1]
        L_dict = {**{j: GO+j for j in range(GO)},
                  **{GO*(i-3)+j: GO*t+GO*(i-4)+j 
                     for j in range(GO)
                     for i in range(4, t+1)},
                  **{GO*(t-2)+j: GO+gs.index(gj*h2) 
                     for j, gj in enumerate(G)},
                  **{GO*(t-2)+GO*(i-3)+j: GO*t+GO*(t-3)+GO*(i-4)+j 
                     for j in range(GO)
                     for i in range(4, t+1)}}
        R_dict = {**{GO*(i-3)+j: GO*(i-1)+j 
                     for j in range(GO)
                     for i in range(3, t+1)},
                  **{GO*(t-2)+GO*i+j: GO*(i+2)+gs.index(gj*hi) 
                     for j, gj in enumerate(G)
                     for i, hi in enumerate(hs[2:])}}
        # This method, exploiting Sage's ability to multiply sparse matrices
        # quickly, is orders of magnitude faster than trying to index and sum. 
        L21 = matrix(BR, nde*GO, ndf*GO, 
                     {(L_dict[j], j):1 for j in range(ndf*GO)}, sparse=True)
        R21 = matrix(BR, nde*GO, ndf*GO, 
                     {(R_dict[j], j):1 for j in range(ndf*GO)}, sparse=True)
        IM = D_cohom*L21*RI_2_diag*(D_cohom*R21).transpose()
    
    # Compute the group representation 
    hom_rep = []
    for g in elms:
        rho_g = matrix(BR, GO, GO, {(gs.index(g*gi), i): 1 
                                    for i, gi in enumerate(G)}, 
                       sparse=True)
        Dg = D_cohom*block_diagonal_matrix(nde*[rho_g.transpose()])
        LD = (Dg*RI_1).change_ring(BR)
        hom_rep.append(LD)

    # Run all consistency checks at the end for code visual simplicity
    # Obviously if the code was stalling at some point one might wich to move
    # these consistency checks earlier. 
    if check_consistency:
        assert G.subgroup(gv).order() == GO
        assert prod(gj for gj in gv) == G.identity()
        assert M1*M2 == 0
        assert M1.right_nullity() - M2.rank() == 2*sigma
        assert D_cohom*N1 == 0
        assert D_cohom.nrows() == 2*sigma
        assert all([LD.is_invertible() for LD in hom_rep])
        if BR(1) == BR(-1):
            print("Warning: if the curve is hyperelliptic the following \
                   tests will produce spurious errors as 1==-1")
        assert all([LD.multiplicative_order()==gi.order() 
                    for LD, gi in zip(hom_rep, elms)])
        assert MatrixGroup(hom_rep).is_isomorphic(G)

        if intersection_matrix:
            assert all([LD*IM*LD.transpose()==IM for LD in hom_rep])

    if intersection_matrix:
        return hom_rep, IM
    else:
        return hom_rep
