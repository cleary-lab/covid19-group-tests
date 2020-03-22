###########################################################################################################
#
# This file contains methods to generate and view group testing designs based on Reed-solomon codes, based on 
# discussions at the IMA group testing workshop, Feb 2012.
#
# --Mary 
#
#####
#
# Some observations:
#
# 1. For any linear code over F_q, if the distance is >= q( 1- 1/r ) then the Kautz-Singleton argument 
#    implies that the corresponding concatenated matrix is r-disjunct. 
#
# 2. For any linear code over F_q, if we mod out by a subspace A with distance \ell, then all cosets
#    have distance \ell as well
#
# 3. In particular, for RS-codes with degree d over F_q, if the order of the matrix is based on cosets 
#    of, say, the leading or constant coefficients,
#    the whole matrix has distance q - d and each sub-block has distance q - d + 1.  
#    Then each block is (d/(d-1))-disjunct and the whole thing was (d+1)/d disjunct.
#    (That is, we have 'modularity').
#
#####
#
#  USAGE:  For example, to get degree-1 Reed-Solomon codes over F_16, ordered lexicographically, do (in the sage shell):
#
#  sage: load rs.sage
#  sage: C,m,polys = create_rs_code(16,1,orderoption="lex", orderfieldelts="lex") 
#  sage: m.show()    // to give you a picture. 
#
#  I can only make 3 pretty pictures with degree 1 codes:
#  1) C,m,polys = create_rs_code(16,1,orderoption="lex", orderfieldelts="lex")  
#            (this is the picture that Or Zuk got)
#            (you actually get the same picture with orderoption="modoutconstant", orderfieldelts="lex", evaluate_backwards=False)
#  2) C,m,polys = create_rs_code(16,1,orderoption="orderorbit", orderfieldelts="mult", evaluate_backwards=False)  
#            (eval_backwards can be T or F for (2), get distinct pretty pictures)
#            This is not strictly modular, but all the stripes go the same way, maybe that's better?
#  3) C,m,polys = create_rs_code(16,1,orderoption="modoutbyleading",orderfieldelts="mult",insideout=True, evaluate_backwards=False) 
#  
#
#
##########################################################################################################


#############################################################################################
# Makes a (concatenated) Reed-Solomon code.
#
# q = order of field
# d = degree of polynomials
# concatenated = True to concatenate with the identity matrix (if false you'll get a matrix over F_q)
# orderoption in {"lex", "modoutconstant","modoutbyleading", "orderorbit",...} see below for descriptions of each ordering scheme
# insideout = True if the concatenation is "inside out," that is, each "level" now corresponds not to a symbol in the original codeword, 
#             but to a row of the binary code (the identity matrix).  This seems like it makes better patterns sometimes, and is
#             just a permutation of the rows
# evaluate_backwards:  If this is false, we map a vector (a_0,...,a_{d-1}) to a polynomial like a_0 + a_1 x + ... + a_n x^n.  If this is
#             true, we do it left-to-right, like a_0 x^n + a_1 x^(n-1) + .....   This has the effect of changing what "lexicographic" means
#
# orderfieldelts is "mult" for a multiplicative ordering or "lex" for a lexicographic ordering.
#
#
# Note that insideout and evaluate_backwards don't add anything over the different orderoptions for linear polynomials, 
# but I think they make a difference for higher order stuff.  In any case, it's nice to explicitly know what's going on. 
#
# Returns: C, picture, polys:
# 	C is the code, either a binary matrix or a matrix of F_q elements, depending on the concatenate option.
#	picture is a plot of the matrix, if it's binary.  To see it, do picture.show().  If C is over F_q, there is no picture.
# 	polys is a list of the polynomials in the order used.  Each polynomial is represented as an element of F_q^(d+1).
######################################################################################################################
def create_rs_code(q,d,concatenated=True,orderoption="lex",insideout=False, evaluate_backwards=True, orderfieldelts="lex",modoutvec=None):
	F.<a> = FiniteField(q)      # make a finite field with generator a of order q
	L = VectorSpaces(F)
	M = L(ZZ^(d+1))             # M =  F_q^(d+1) is where the polynomials live 
	## Note: another option is to go in the polynomial ring instead of F_q^(d+1):	
	## R.<t> = PolynomialRing(F) gives the ring of polynomials over F, in indeterminant t
	## R.polynomials(max_degree=d) gives a list of all the polynomials of appropriate degree
	## but I think I want to use sage's vector space structure without coercing each time
	# fix an ordering on the elements
	if order(F).is_prime():
		elts = [a*s for s in range(q)]  # sage gives the generator a=1 if the additive group is cyclic
	else:
		if orderfieldelts=="lex":
			elts = order_field_elements(F)
		elif orderfieldelts=="mult":
			elts = [F.zero_element()] + [ a^s for s in range(1,q) ]  # else a is the multiplicative generator
		else:
			print("Bad option for orderfieldelts (must be 'lex' or 'mult')!")
			exit(0)
	## lots of different ways to generate ordering of polynomials:
	if orderoption is "lex":
		## 1: Lexicographic according to my previous ordering
		polys = get_lexicographic( elts, M )
	elif orderoption is "modoutconstant":
		## 2: mod out by the subspace generated by the constant term, lexicographically in there
		polys = order_polys_mod_out_subspace( M, M.subspace( [M.basis()[0]] ),elts)  
	elif orderoption is "modoutbyleading":
		## 3: mod out by the subspace generated by the leading term, lexicographically within there
		## this is what STD does when q is prime.
		polys = order_polys_mod_out_subspace( M, M.subspace( [M.basis()[-1]] ),elts)  
	elif orderoption is "orderorbit":
		## 4:  This fails to be modular, but it has really nice diagonals.  
		## This is taking all the monic polynomials and creating chunks out of their orbits under mult by a generator of F_q^*.
		polys = order_by_orbit( M )
	elif orderoption is "modout":
		# modoutvec should be a vector of indices into elts.
		if modoutvec == None:
			print("I need a vector")
			exit(0)
		c = copy(M.basis()[0])
		for i in range(len(modoutvec)):
			c.set(i,elts[modoutvec[i]])
		print("I'm about to mod out by", c)
		polys = order_polys_mod_out_subspace( M, M.subspace( [c] ),elts )
	else: 
		print("Bad order option")
		exit(0)
	M = matrix( [[ apply_poly(p, elt,evaluate_backwards) for p in polys] for elt in elts ] ) # make the (non-concatenated) RS code
	if not concatenated:
		return M, None, polys
	C = concatenate(M,q,polys,elts,insideout)
	# for printing purposes, organize C by qxq blocks (or (q-1)xq blocks for the orbit scheme).
	if orderoption is "orderorbit":
		C.subdivide([q*i for i in range(C.nrows()/q)], [1] + [ 1 + (q-1)*i for i in range( (C.ncols()-1)/(q-1)) ] )
	else:
		C.subdivide( [ q*i for i in range(C.nrows()/q) ], [q*i for i in range(C.ncols()/q) ] )
	picture = matrix_plot(C, subdivisions=True, subdivision_style=[ dict(color="orange",thickness=.5), dict(color="red", thickness=.5) ],fontsize=7  )
	if q >= 16:
		picture.SHOW_OPTIONS['dpi'] = 600 # it looks gross with the default dpi=100
	return C, picture ,polys  # to plot, do picture.show() 


# apply a polynomial p in F_q^d to x in F_q
def apply_poly( p, x , evaluate_backwards=False):
	# if we have (a0 a1 a2 a3 ... ad) then we evaluate a0 + a1x + ... + ad x^d. 
	if not evaluate_backwards:
		return p[0] + sum([ p[s] * x^s for s in range(1,len(p)) ])
	else:
		# now do ad + a(d-1) x + ...  This is the same as doing the lex order backwards.
		d = len(p)
		return p[d-1] + sum([ p[d-s-1] * x^s for s in range(1,d) ])

# concatenates an RS code with the identity code to get a binary matrix
def concatenate(M,q, polys, elts, insideout):
	C = matrix( [ [ 0 for i in range (M.ncols()) ] for j in range(q^2) ] )
	for i in range(len(elts)):  # go through the layers i=1,...,q
		for j in range(len(elts)): # the place within the layer
			for k in range( M.ncols()) : # go across the matrix
				if M[i,k] == elts[j]:
					if not insideout:
						C[q*i + j,k] = 1
					else:
						C[i + q*j, k]=1
	return C

################################################
# Order the field elements lexicographically
################################################
def order_field_elements(F):
	ret = []
	current = F.zero_element()
	order_field_helper(F, ret, current, F.polynomial().degree()-1)
	return ret

def order_field_helper(F, ret, current, index):
	a = F.gen()
	if index < 0:
		ret.append(copy(current))
		return
	for b in range(F.characteristic()):
		current += b*a^index
		order_field_helper( F, ret, current, index - 1)
		current -= b*a^index

############################################
# Methods for ordering polynomials
############################################

# given an ordering on the elements, order the elements of M = F_q^(d+1) lexicographically
def get_lexicographic( elts, M ):
	ret = []
	current = [ None for i in range(M.dimension()) ]
	get_lex_helper(elts, M, ret, current, 0 )
	return ret

def get_lex_helper(elts, M, ret, current, index):
	if index == M.dimension():
		u = copy(M.zero_element())
		for i in range(len(current)):
			u.set(i, current[i] )
		ret.append(u)
		return
	for e in elts:
		current[index] = e
		get_lex_helper( elts, M, ret, current, index + 1)

# A is a subspace of F_q^(d+1) of order q, return the cosets (lexicographic ordering within cosets)
def order_polys_mod_out_subspace( M, A, elts, submod=None ):
	Q = M.quotient(A)
	# each element of Q gets its own block.  But we need to figure out how to order the block.
	ret = []
	liftmap = Q.lift_map()
	# now we need to order the elements of Q.  Either go lexicographically or recurse and use the passed-in subspace.
	if submod is None:
		if dimension(Q) > 1:
			real_order = A.list()
		else:
			# order the elements of A to line up with the elements of elts.
			lex_order = get_lexicographic(elts, M)
			tmporder = []
			for y in A.list():
				j = lex_order.index( y )
				tmporder.append( (j,y) )
			tmporder.sort()
			real_order = []
			for x in tmporder:
				real_order.append(x[1])
	else:
		real_order = order_polys( Q, submod, None ) # right now I only support one level of recursion...
	for coset in Q.list():
		for elt in real_order:
			ret.append( liftmap(coset) + elt )
	return ret

# each chunck corresponds to a monic polynomial, and the (q-1) elements in each chunk are the orbit of that poly under mult by a.
def order_by_orbit( M ):
	a = M.base_field().gen()
	q = M.base_field().order()
	A = M.subspace( [M.basis()[-1]] ) # mod out by the leading term
	Q = M.quotient(A)
	ret = [ M.zero_element() ]
	# first do the constants
	for coset in Q.list():
		if coset == Q.zero_element():
			continue
		ret.append( Q.lift(coset) )
	# now do all of the monic polys 
	for coset in Q.list():
		if coset == Q.zero_element():
			continue
		monic_poly = M.basis()[-1] + Q.lift(coset)
		for s in range(1, q):
			ret.append( a^s * monic_poly )
	return ret

##################################################
# Methods for checking distance (sanity checks)
#################################################

# get the distance of the code C, restricted to cols if not None
def get_dist( C, cols=None ):
	if cols is not None:
		index = cols
	else:
		index = range(C.ncols())
	dist = C.nrows()
	for i in range(len(index)):
		for j in range(i+1, len(index)):
			diff = 0
			for k in range(C.nrows()):
				if C[k,index[j]] != C[k,index[i]]:
					diff += 1
			if diff < dist:
				dist = diff
	return dist

# make sure that each contiguous q-sized sub-block of C has good distance
def sanity_check_dist( C, q ):
	total_dist = get_dist(C)
	worst_sub_dist = C.nrows()
	for i in range(C.ncols()/q):
		curcols = range(q*i,q*(i+1))
		sub_dist = get_dist(C, cols=curcols)
		if sub_dist < worst_sub_dist:
			worst_sub_dist = sub_dist
	return total_dist, worst_sub_dist
	
#####################################################
# code for playing around....
####################################################

## are there any subspaces for the F_16 linear case that play nice with multiplicative structure?
def test_subspaces(insideout=False, evaluate_backwards=True):
	pictures = []
	for v in [[1,0],[0,1]]:
		C, m, polys = create_rs_code( 16, 1, orderoption="modout", insideout=insideout, evaluate_backwards=evaluate_backwards, orderfieldelts="lex", modoutvec=v )	
		pictures.append(m)
	return pictures
	