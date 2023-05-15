### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ 93e62420-39d1-4813-9bd5-7fa12286d190
begin
	using StructArrays 	# For array of corner structs
	using BenchmarkTools # For benchmarking
	using Random 		# For randomly generating points
end

# ╔═╡ 3138335c-5ce5-4f70-bf3e-7d093ce8de13
md"# Delaunay Triangulation in 3D"

# ╔═╡ 9f95fab7-f3d3-46c1-ac32-537b84196a35
md"### Settings and Global Utilities"

# ╔═╡ 23f069f2-c3a8-47ff-8af5-1f3a242e6bb4
# Constants
begin
	const SIMPperV::Int32 = 100 	# max simplices per vert
	# const FORGETLIMIT = 10000 	# forget simplicesafter this many (for better locality of ref? Doesn't help  - from delaunay3.h)

	const STACKMAX::Int32 = 2000

	const HIGHCOORD::Int32 = 0x400 #0x4000
	const COORDMASK::Int32 = 0x3FF #0x3FFF

	# const CornerIndex = Int32 	# TODO: Use this for clarity when dealing with corner indices?
end;

# ╔═╡ 512b21a0-1e6f-44b8-a05c-38d45ec7467f
md"### Points and Corners"

# ╔═╡ cada11d4-0aa9-4b74-a2ff-593b2497e9f6
#= 
A Point contains (x,y,z) spacial coordinates and a lifted sq coord. 
=#
struct Point
	#index::Int32 	# Index back to original line in file
	x::Int32
	y::Int32
	z::Int32
	sq::Int64
	
	Point(x, y, z) = new(x, y, z, x*x+y*y+z*z)
	Point(x, y, z, sq) = new(x, y, z, sq)
end

# ╔═╡ 083f724b-2009-4a49-ae1e-61363c3d7165
#=
A Corner contians vertex location (Point) as well as the index of the opposite corner in the neighboring simplex
=#
mutable struct Corner
	vertex::Point
	opp::Int32  	# Index of opp corner in corner vector
	
	Corner(v, opp) = new(v, opp)
	Corner(v) = new(v)	
	Corner() = new()
end

# ╔═╡ 2e9a743d-78d2-4a89-a400-3f6ade10b7e0
# Utility funcitons for handling points or corners
begin
	# TODO: Should this also check sq? (In case we have the origin? It would match with point @ infinity)
	@inline function equalV(a::Point, b::Point)::Bool
		return a.x == b.x && a.y == b.y && a.z == b.z
	end
	
	# Is the given point the point @ infinity?
	@inline function infiniteV(pv::Point, vert::Point)::Bool
		return pv == vert
	end
end;

# ╔═╡ c368ce56-b068-45ac-979f-084e3175c8e0
md"## Using Array of Corners"

# ╔═╡ 75695a4f-7c0c-4879-88d8-c2fc7acbccca
md"Every 4 corners make up one simplex. This makes it easy to find the simplex a corner is in. The 'this.corners' array later will leave the first 3 indices unused. This way, simplex #1 starts at corner #4, and the starting corner of a simplex is that simplex's number times 4."

# ╔═╡ 81ce6604-91ef-4751-b4ff-1535dad19c40
# Simpex manipulation tables
begin
	# c+offset[i][INDEX(c)+1] advances c to (c+i)mod5
	const offset::Array{Int8, 2} = [
		0 0 0 0;
		1 1 1 -3;
		2 2 -2 -2;
		3 -1 -1 -1 
	]

	# drop[i] contains new vertex order after vertex i is dropped and replaced by pv on same side
	# const drop::Array{Int8, 2} = [
	# 	2 1 3;
	# 	0 2 3;
	# 	1 0 3;
	# 	0 1 2
	# ]
	const drop::Array{Int8, 2} = [
		3 2 4;
		1 3 4;
		2 1 4;
		1 2 3
	]

	# offdr[i] contains drop(i) - INDEX(i) - 1
	# const offdr::Array{Int8, 2} = [
	# 	2 1 3;
	# 	-1 1 2;
	# 	-1 -1 1;
	# 	-3 -2 -1
	# ]
	const offdr::Array{Int8, 2} = [
		 2  1  3;
		-1  1  2;
		-1 -1  1;
		-3 -2 -1
	]

	# invdrop[i][k] = j whenever drop[i][j] = k. 
	# 4's signal i=k; bad because i is dropped
	# const invdrop::Array{Int8, 2} = [
	# 	4 1 0 2;
	# 	0 4 1 2;
	# 	1 0 4 2;
	# 	0 1 2 4
	# ]
	# Added 1 to handle Julia's 1-indexing
	const invdrop::Array{Int8, 2} = [
		5 2 1 3;
		1 5 2 3;
		2 1 5 3;
		1 2 3 5
	]
end;

# ╔═╡ 52f11677-636c-4a90-9eeb-cf5e4b948bfe
# Macros taken from d3.h
# Note: Julia uses 1-indexed arrays
# Corner indeces 1, 2, 3 unused. Simplex 1 starts with corner 4. Corners are 0-indexed in a simplex

begin
	# Calculate a mod 4
	@inline function MOD4(a)::Int32
		return a & 0x3
	end

	# Find the simplex ID of the given corner
	@inline function SIMPLEX(cornerIdx::Int32)::Int32
		return cornerIdx >> 2
	end

	# What is the index of the given corner in its simplex?
	# Indices range from 0:3
	@inline function INDEX(cornerIdx::Int32)::Int32
		return MOD4(cornerIdx)
	end
	
	# Get the specified corner of the given simplex
	@inline function CORNER(simp::Int32, index::Int32)::Int32
		return (simp << 2) + index
	end
	# With different int types
	@inline function CORNER(simp::Int64, index::Int64)::Int32
		return CORNER(Int32(simp), Int32(index))
	end

	# Find the base corner (first corner) of the simplex containing the given corner
	@inline function BASECORNER(cornerIdx::Int32)::Int32
		return (cornerIdx & 0xFFFFFFFC)
	end

	# Find the last corner of the simplex containing the given corner
	@inline function LASTCORNER(cornerIdx::Int32)::Int32
		return (cornerIdx | 0x3)
	end

	# Cycles through the [ 1, 1, 1, -3 ] array to iterate through corners of a simplex
	@inline function INCREMENT(corner::Int32)::Int32
		return @inbounds corner + offset[2, INDEX(corner)+1]
	end

	# Generates one random bit (0 or 1)
	@inline function RANDBIT()::Int8
		return rand((0,1))
	end
end;

# ╔═╡ f92bf3ea-b573-4790-8b84-348ddef64aca
md"# Representing Spheres"

# ╔═╡ 0344f8c1-49d9-482b-b9e0-60ca3bfb8bab
# Save sphere questions
struct Sphere
	x::Int64
	y::Int64
	z::Int64
	sq::Int64
end

# ╔═╡ 8ddc46f4-2e4b-4184-9a49-5c5d6d3be112
# Sphere utility functions
begin
	# 2x2 determinant of [[a, b], [c, d]]
	@inline function det2(a, b, c, d)::Int64 	# Is this the right type to use? does it matter?
		return a * d - b * c
	end

	# Instead of storing sp.w, this can be calculated using another point (sv)
	@inline function spdot(sp, pv, sv)
		sp.x*(pv.x-sv.x) + sp.y*(pv.y-sv.y) + sp.z*(pv.z-sv.z) + sp.sq*(pv.sq-sv.sq)
	end
end;

# ╔═╡ bf22c408-6d8d-45d3-a4d7-cbda3dd9fa1a
md"### makeSphere and inSphere"

# ╔═╡ 8ced7d59-909a-46b1-954e-0c91324ec2d1
#= Generate a sphere from 4 given points
		pv is 'point of interest'
		vert represents the first element in the verts array (point @ infinity)
=#
# Note: d3.c code was void and took a sphere pointer as input. Here, it is a return value
function makeSphere(v0::Point, v1::Point, v2::Point, pv::Point, vert::Point)::Sphere
	if (!infiniteV(v0, vert))
		x0 = v0.x - pv.x; y0 = v0.y - pv.y;
		z0 = v0.z - pv.z; sq0 = v0.sq - pv.sq;
	else 
		x0 = v0.x; y0 = v0.y;
		z0 = v0.z; sq0 = v0.sq;
	end

	if (!infiniteV(v1, vert))
		x1 = v1.x - pv.x; y1 = v1.y - pv.y;
		z1 = v1.z - pv.z; sq1 = v1.sq - pv.sq;
	else 
		x1 = v1.x; y1 = v1.y;
		z1 = v1.z; sq1 = v1.sq;
	end

	x2 = v2.x - pv.x; y2 = v2.y - pv.y;
	z2 = v2.z - pv.z; sq2 = v2.sq - pv.sq;

	# 2x2 minors
	xy = det2(x0, x1, y0, y1)
	xz = det2(x0, x1, z0, z1)
	yz = det2(y0, y1, z0, z1)
	xs = det2(x0, x1, sq0, sq1)
	ys = det2(y0, y1, sq0, sq1)
	zs = det2(z0, z1, sq0, sq1)

	x = -y2*zs + z2*ys - sq2*yz
	y =  x2*zs - z2*xs + sq2*xz
	z = -x2*ys + y2*xs - sq2*xy
	sq = x2*yz - y2*xz + z2*xy

	# Store 2x2 minors
	sp = Sphere(x, y, z, sq)

	return sp
end

# ╔═╡ 9db9a542-ed11-4408-bbe9-a2b1d56bca96
# Is the point pv inside sphere sp?
@inline function inSphere(sp::Sphere, pv::Point, sv::Point)::Int64
	d::Int64 = spdot(sp, pv, sv) 	# Inside if negative
	
	return d
end

# ╔═╡ bd61ed12-f3c8-452d-b848-b4d9994164af
md"# Stack"

# ╔═╡ c21855ca-c54d-4e4b-b68c-b3b548eda1db
# basic stack
begin 
	const StackEltType = Int32 # stack entry type
	
	mutable struct Stack
		const st::Vector{StackEltType} 	# stack space
		sp::Int32 					# stack pointer
		Stack() = new(zeros(StackEltType,STACKMAX), Int32(0))
	end

	@inline	function pop!(st::Stack)::StackEltType
		@inbounds local ret = st.st[st.sp] #pop
		st.sp-=oneunit(st.sp)
		ret
	end

	@inline	push!(st::Stack, x::StackEltType) = @inbounds 	st.st[st.sp+=oneunit(st.sp)] = x # push

end;

# ╔═╡ d52396af-0dfb-4fa0-9942-fba52d282eae
# Utility functions for the stack
begin
	@inline function isEmpty(st::Stack)::Bool
		return st.sp < 1
	end
end;

# ╔═╡ 27daafcc-f7b1-43bc-8689-70fef984fa5f
md"### Testing Stack"

# ╔═╡ 3f8be73a-4ce9-4ef1-bb94-e572e82ac382
begin
	function testStack(xlist::Vector{StackEltType})
		st = Stack()
		
		for x in xlist
			push!(st, x)
		end
		
		for x in (xlist)
			pop!(st) #@assert(pop!(st) == x)
		end
		(st.sp,st.st)
	end
	
	# testStack(rand(StackEltType,10))
	@benchmark testStack(xlist) setup=(xlist=rand(StackEltType, STACKMAX))
end

# ╔═╡ 08e8aeb5-8d8e-44b7-bdd3-1f8e036611b6
md"# Initialize"

# ╔═╡ f4fff776-1f97-4bd9-9fd5-7d22780a7028
# Struct to hold state information
# Generally referred to as 'this'
mutable struct StateType
	verts::Vector{Point} 	# vertices; 0th is point at infinity!
	corners::StructVector{Corner} # corners (called 's' in delaunay3.h)
	sphs::Vector{Sphere} 	# spheres

	active::Vector{Int32} 	# which spheres are active 
							# flag: -1 unused, 0 dead, 1 alive
							# TODO: Should this be a different type? (like an enum)
							# TODO: Why not make this an attribute of the sphere?
	
	freeSimp::Int32		# head for free list for simplices 
	liveSimp::Int32 	# latest tetra; known to be live
	maxSimps::Int32 		# AUDIT Only. The number of simplices (active or not)
	limitmaxSimps::Int64 	# limit on # of created simplices, spheres, & corners/4
	
	dfs::Stack				# DFS stack to find dead simplices
	idfs::Stack				# DFS stack for simplices adj to infinite vertex
	nhbr::Stack				# stack for dead corners with live neighbors
	kill::Stack				# stack for base corners of simplices to recycle

	# Constructor
	StateType() = new()
end

# ╔═╡ 77719559-7da6-49ac-b5f3-e5de2428bb8a
# TODO: Fix this comment. Indices are different here.

# The table indoff has the index offset for where 
# the vertex of index i will be found in the tetrahedron opposite c.  
# Vertex s[BASECORNER(c)+i].v will be at s[s[c].opp + indoff[INDEX(c)][INDEX(s[c].opp)][i]].v 
# (except that s[c].v and s[s[c].opp]].v are different, and opposite sides of common pl).
# Note that indoff[*][j] uses each offset -j:4-j exactly once, 
# that indoff[i][j][i] = 0 for all possible i and j (so s[c].v is opposite of s[s[c].opp].v), 
# and that the tetra will never change, TETRA(s[c].opp) == TETRA(CORNERINOPP(i,c)).


# This is the table of where the indices go; used in ASSERTS only.
# Replacing  \ With V falling at position c.opp:
#   corner c  \0ABCD 1ABCD 2ABCD 3ABCD 
#       A    0: xxxx  BVCD  CBVD  BCDV
#       B    1: VACD  xxxx  ACVD  CADV 
#       C    2: vBAD  AVBD  BAVD  ABDV 
#       D    3: Vabc  bvac  ABVC  BACV
#
# Offsets from c.opp  I=-1, Z=-2, B = -3, H = -4
# Replacing \0ABCD 1ABCD 2ACBD 3ABCD 
#     A    0: xxxx  0I12  0IZ1  0BZI 
#     B    1: 1023  xxxx  Z0I1  Z0BI 
#     C    2: 1203  I102  IZ01  BZ0I 
#     D    3: 1230  1I20  ZI10  ZBI0 
# Get these by subtracting 01234 from columns in order
begin
	function CORNERINOPP(this::StateType, i::Int64, c::Int32)
		return this.corners.opp[c] + indoff[INDEX(this.corners.opp[c])+1, i+1, INDEX(c)+1]
	end
	
	const indoff::Array{Int8, 3} = cat(
		[ 5  5  5  5;  
		  0 -1  1  2;  
		  0 -1 -2  1;  
		  0 -3 -2 -1],
		
		[ 1  0  2  3;  
		  5  5  5  5; 
		 -2  0 -1  1; 
		 -2  0 -3 -1],
		
		[ 2  1  0  3; 
		 -1  1  0  2; 
		 -1 -2  0  1; 
		 -3 -2  0 -1],
		
		[ 1  2  3  0;  
		  1 -1  2  0; 
		 -2 -1  1  0; 
		 -2 -3 -1  0],
		
		dims=3
	)
end;

# ╔═╡ 25fb3749-b637-4095-b277-8b3325161a25
# Utility functions for dealing with StateType 'this'
begin
	# Kill simplex p
	@inline function KILL(this::StateType, p::Int32)
		@inbounds this.active[p] = -1
	end

	# Is this a dead or killed simplex?
	@inline function DEAD(this::StateType, p::Int32)
		@inbounds return this.active[p] <= 0
	end

	# Allocate or free space for simplex
	# When allocating, this.liveSimp is the location of the new simplex
	@inline function STARTSIMP(this::StateType)
		# startTetraCnt+=1
		this.liveSimp = this.freeSimp

		if (0 > this.corners.opp[CORNER(this.liveSimp, Int32(0))])
			println("MaxSimps: ", this.limitmaxSimps)
			println("FreeSimp: ", this.freeSimp)
			println(CORNER(this.liveSimp, Int32(0)))
			@inbounds println(this.corners[CORNER(this.liveSimp, Int32(0))])
			println()
		end
		
		@inbounds this.freeSimp = this.corners.opp[CORNER(this.liveSimp, Int32(0))]
		@assert DEAD(this, this.liveSimp) "Reusing existing simplex?"
		@inbounds this.active[this.liveSimp] = 1

		if (this.maxSimps < this.liveSimp) 
			this.maxSimps = this.liveSimp
		end
		# AUDIT
		if (this.maxSimps >= this.limitmaxSimps)
			print("AUDIT: ")
			print(this.liveSimp)
			println(" > limitmaxSimps")
		end
	end

	@inline function FREESIMP(this::StateType, p::Int32)
		@inbounds @assert this.active[p]<0 "Freeing already free simplex?"
		@inbounds this.active[p] = 0
		@inbounds this.corners.opp[CORNER(p, Int32(0))] = this.freeSimp
		this.freeSimp = p

		# If using FORGETLIMIT
		# if ( (this.maxSimps - p) < FORGETLIMIT )
		# 	this.corners.opp[CORNER(p, Int32(0))] = this.freeSimp
		# 	this.freeSimp = p
		# end
	end

	# Set corner's vertex and opposite in simplex structure
	@inline function setCornerVC!(this::StateType, c::Int32, vertex::Point, opp::Int32)
		@inbounds this.corners.vertex[c] = vertex
		@inbounds this.corners.opp[c] = opp
	end
	# Also sets the opposite's opp to this corner
	@inline function setCornerVCN!(this::StateType, c::Int32, vertex::Point, opp::Int32)
		setCornerVC!(this, c, vertex, opp)
		@inbounds this.corners.opp[opp] = c
	end
end;

# ╔═╡ d7cdc261-232c-4560-95f1-52aa2e4074c7
md"initialOpp array:"

# ╔═╡ 894587ef-8eee-4eeb-8354-2dbbda06a852
# When initializing, this is what gets filled in
	# recall CORNER() takes a simplex number (>=1) and 0<= index <=3 in that simplex and returns the index of that corner in the array of corners

#=
C code: 
static const int initialopp[5][4] = {
  {CORNER(1,1), CORNER(2,0), CORNER(3,1), CORNER(4,0)}, 
  {CORNER(2,1), CORNER(0,0), CORNER(3,0), CORNER(4,1)},  
  {CORNER(0,1), CORNER(1,0), CORNER(3,2), CORNER(4,2)}, 
  {CORNER(1,2), CORNER(0,2), CORNER(2,2), CORNER(4,3)}, 
  {CORNER(0,3), CORNER(1,3), CORNER(2,3), CORNER(3,3)}};
=#
# Added 1 to the simplex numbers since they are 1-indexed

const initialOpp::Array{Int32,2} = [
	CORNER(2,1) CORNER(3,0) CORNER(4,1) CORNER(5,0); 
	CORNER(3,1) CORNER(1,0) CORNER(4,0) CORNER(5,1);
	CORNER(1,1) CORNER(2,0) CORNER(4,2) CORNER(5,2);
	CORNER(2,2) CORNER(1,2) CORNER(3,2) CORNER(5,3);
	CORNER(1,3) CORNER(2,3) CORNER(3,3) CORNER(4,3)
];

# ╔═╡ 492d064a-fc8e-45d6-a548-e206d7a6e65b
md"initialV array:"

# ╔═╡ b112c2d3-a38e-41e0-a2be-839faa770b8e
#=
C code:
static const int initialv[2][5][4] = {{{1,2,3,4}, {2,0,3,4}, {0,1,3,4}, {1,0,2,4}, {0,1,2,3}},
                                      {{0,2,3,4}, {2,1,3,4}, {1,0,3,4}, {0,1,2,4}, {1,0,2,3}}};
=#
# Added 4 to these b/c of 1-indexing in corners array and skipping 1, 2, and 3

# TODO: Verify these

const initialV::Array{Int32, 3} = cat([
	2 3 4 5;
	3 1 4 5;
	1 2 4 5;
	2 1 3 5;
	1 2 3 4
],
[
	1 3 4 5;
	3 2 4 5;
	2 1 4 5;
	1 2 3 5;
	2 1 3 4
], dims=3);

# ╔═╡ 9c700c27-565f-4a9b-b031-229f434d0512
md"### Initialize!"

# ╔═╡ 5313f41e-19d7-4257-82e0-d04db110050b
#= 
Initialize the state from the given array of points. Allocates memory for corners and spheres, sets up free list for simplices, and creates spheres for the first five points.
REQUIRES: vertArray[1] contains the point @ infinity and vertArray[2:5] contain four finite points; These first five points must be in general posiiton
=#

function d3Initialize!(this::StateType, vertArray::Vector{Point}, nVert::Int32)::Bool
	# init state
	this.verts = vertArray
	this.limitmaxSimps = SIMPperV * nVert 	# space for spheres and corners
	
	this.sphs = Vector{Sphere}(undef, this.limitmaxSimps)
	this.active = Vector{Int32}(undef, this.limitmaxSimps)

	# For some reason adding this temp 'corners' variable fixed a bug...
	corners = StructVector{Corner}(undef, 4 * this.limitmaxSimps)
	this.corners = corners;

	# init simplices
	last::Int64 = -1
	this.freeSimp = this.limitmaxSimps

	# Set up free list of tetrahedra
	while (this.freeSimp > 6)
		this.freeSimp -= 1
		@inbounds this.active[this.freeSimp] = 0 	# starts dead
		
		@inbounds corners.opp[CORNER(this.freeSimp, zero(Int32))] = last
		last = this.freeSimp
	end

	# Make first sphere
	# TODO: Change these indices to be more clear, they get overridden anyway
	@inbounds this.active[5] = 2 	#TODO: Why is this 2?
	@inbounds this.sphs[5] = makeSphere(this.verts[1], this.verts[2], this.verts[3], this.verts[4], this.verts[1])	# this.verts[1] is point @ infinity
	@inbounds d = spdot(this.sphs[5], this.verts[5], this.verts[4]) 	# if d<0, then we need to swap (i.e. use the second InitialV lookup table)
	
	if (d == 0)
		# First five points not in general position
		println("ERROR: Need first five vertices to be in general position")
		return false
	end

	# Make sure this 'infinite sphere' has the proper orientation
	#= Pretty sure this code from d3.c is useless (sphs[5] gets overwritten)
	if (d < 0)
		# this.sphs[5].x = -this.sphs[5].x
		# this.sphs[5].y = -this.sphs[5].y
		# this.sphs[5].z = -this.sphs[5].z
		# this.sphs[5].sq = -this.sphs[5].sq
		this.sphs[5] = Sphere(-this.sphs[5].x, -this.sphs[5].y, -this.sphs[5].z, -this.sphs[5].sq)
	end
	=#

	# Create the first 5 simplices
	for p in 1:5
		# Make the 4 corners
		for j in 0:3
			# Pay attention to orientation when assigning vertices
			@inbounds nextVert::Int32 = initialV[p, j+1, (d<0)+1]
			@inbounds opp::Int32 = initialOpp[p, j+1]
			@inbounds corners[CORNER(p, j)] = Corner(this.verts[nextVert], opp)
		end

		@inbounds this.sphs[p] = makeSphere(corners.vertex[CORNER(p, 0)], 
								   corners.vertex[CORNER(p, 1)],
								   corners.vertex[CORNER(p, 2)],
								   corners.vertex[CORNER(p, 3)],
								   this.verts[1])
		
		@inbounds this.active[p] = 1
	end
	
	this.liveSimp = 5 
	this.maxSimps = 5

	return true
end

# ╔═╡ f5b545bc-ab45-4c6c-9843-c076b79d4915
md"# Locate Sphere"

# ╔═╡ f3d6901e-8777-4505-91c6-4b9be755f5c3
#=
Location routine which walks a mesh stored in corners, active, and sph starting from SIMPLEX start to find point pv.
Result is the index of the simplex whose sphere strictly contains pv.
(The simplex itself need not contain pv.)
Return code is positive if we succeed (# of location steps at present)
Return code is 0 if we failed, perhaps because points with small weight have no Voronoi cells
Return code is -1 if we found the duplicate of a vertex in the mesh. IN THIS CASE, result is the location of the corner where we found the duplicate
=#

# Note: d3.c code took a result pointer as input. Here, it is part of the return value
function d3LocSphere(this::StateType, pv::Point, start::Int32)::Tuple{Int32, Int32}

	result::Int32 = 0
	guard::Int32 = 2 * this.maxSimps + 4 	# prevent infinite loops

	c1::Int32 = CORNER(start, Int32(0)) 	# corner in start
	@inbounds s1::Sphere = this.sphs[start] 			# sphere at start
	@inbounds I1::Int64 = inSphere(s1, pv, this.corners.vertex[LASTCORNER(c1)]) 	# check if inside start sphere

	if (I1 < 0) # Found one on first try!
		result = start
		return (start, 1)
	end
	
	while ((guard -= 1) > 0) 

		@inbounds c2::Int32 = this.corners.opp[c1]
		@inbounds s2::Sphere = this.sphs[SIMPLEX(c2)]
		@inbounds I2::Int64 = inSphere(s2, pv, this.corners.vertex[LASTCORNER(c2)])

		if (I2 < 0) # Found one!
			result = SIMPLEX(c2)
			steps::Int32 = 2*this.maxSimps + 5 - guard  	# number of steps
			return (result, steps)
		end
		
		d = s2.sq * I1 - s1.sq * I2 	# Warning: if s1 & s2 are same sphere, this is zero
		# locateSideCnt += 1
		if (d == 0) # may be on two spheres -- check for duplicate vertex
			if (@inbounds equalV(pv, this.corners.vertex[c2]))
				 result = c2
				# AUDITDROP # TODO: ?
				return (result, -1)
			end

			# Use offset lookup table to check each other vertex in this simplex
			j = INDEX(c2) + 1
			
			if (@inbounds equalV(pv, this.corners.vertex[c2+offset[1,j]]))
				@inbounds result = c2 + offset[1,j]
				return (result, -1)
			end
			if (@inbounds equalV(pv, this.corners.vertex[c2+offset[2,j]]))
				@inbounds result = c2 + offset[2,j]
				return (result, -1)
			end
			if (@inbounds equalV(pv, this.corners.vertex[c2+offset[3,j]])) 
				@inbounds result = c2 + offset[3,j]
				return (result, -1)
			end
			
			# otherwise, no duplicate; we probably have s1 == s2. (Rare in protein data)
			# Choose random direction, as a hack. (Should do a plane computation)
			d = rand((-1, 1))
		end

		if (d < 0) # on I1 side
			c1 = INCREMENT(c1)
		else # on I2 side
			c1 = INCREMENT(c2)
			s1 = s2
			I1 = I2
		end
	end

	result = 0
	return (result, 0)
end

# ╔═╡ 9f387d99-7848-4bbb-83e3-160f6fe08dfb
md"# Inserting a vertex into the structure"

# ╔═╡ f9b13d44-ecde-4133-86de-bdb9506dd99a
md"### DFS!"

# ╔═╡ 4de76870-e6c6-48de-a954-75e751ef2944
# DFS
@inline function DFS!(this::StateType, pv::Point, p::Int32)
	# DFS
	while (!isEmpty(this.dfs))
		c::Int32 = pop!(this.dfs)
		s = SIMPLEX(c)

		@inbounds @assert DEAD(this, SIMPLEX(this.corners.opp[c])) "dfs stack element with non-dead neighbor"

		if (DEAD(this, s)) 	# dead already
			continue
		end

		# Is pv in, out, or on pth sphere?
		@inbounds d::Int64 = inSphere(this.sphs[s], pv, this.corners.vertex[LASTCORNER(c)])
		if (d < 0)
			# println("in")
			# Kill and continue DFS if pv is strictly inside
			KILL(this, s)
			push!(this.kill, s)

			# Stack neighbors to check
			j::Int32 = INDEX(c) + 1
			@inbounds push!(this.dfs, this.corners.opp[c + offset[2, j]])
			@inbounds push!(this.dfs, this.corners.opp[c + offset[3, j]])
			@inbounds push!(this.dfs, this.corners.opp[c + offset[4, j]])
			
		elseif (d > 0 || this.sphs[p].sq > 0)
			# println("out")
			# pv is outside (or on with sp finite), 
			# so c is live neighbor of dead opp this.corners[c].opp
			@inbounds push!(this.nhbr, this.corners.opp[c])	# remember old corner to hook simplex into mesh later
			STARTSIMP(this)		# make new simplex liveSimp
			newb = CORNER(this.liveSimp, Int32(3))	# last corner of new simp
			setCornerVCN!(this, newb, pv, c) 	# last corner is pv; also set opposite corner c. Do rest later (TODO: Do we do the rest later?)

		else 
			# println("on")
			# d == 0 && this.sphs[p] is infinite: handle two special cases
			if (@inbounds this.sphs[SIMPLEX(this.corners.opp[c])].sq == 0)
				# Then if c stays alive, we make simp to it (flat, but infinite)
				push!(this.idfs, c)
				
			else
				# dead sphere is finite; kill c and make simps to neighbors, if they stay alive
				KILL(this, p) 			# kill and
				push!(this.kill, p)	# stack simp

				# Stack neighbors to check
				j = INDEX(c) + 1
				@inbounds push!(this.idfs, this.corners.opp[c + offset[2, j]])
				@inbounds push!(this.idfs, this.corners.opp[c + offset[3, j]])
				@inbounds push!(this.idfs, this.corners.opp[c + offset[4, j]])
				
			end
		end
	end
end

# ╔═╡ fcf8c6cd-d50d-4991-a5dd-5fc0bcc8483b
md"### newSimpCorner!"

# ╔═╡ dd5593da-8154-4988-820f-a73a24c9a3d6
function newSimpCorner!(this::StateType, pv::Point, newb::Int32, dead::Int32, vNum)::Point
	jdead::Int32 = INDEX(dead) 		# index of dropped corner
	
	# TODO: Maybe fix the drop table values to avoid the unnecessary -1
	# TODO: Also fix indices that have +1's 
	j::Int32 = jdead
	@inbounds i::Int32 = drop[j+1, vNum]-1 	# old index of new corner 1
	c::Int32 = dead + i 		# Note: i = INDEX(c)
	@inbounds v::Point = this.corners.vertex[c] 	# copy vertex v0
	@inbounds nc::Int32 = this.corners.opp[c] 		# go to neighbor

	# In simplex opp c, find new location of j. That is new c. New j = INDEX(c.opp) + 1.
	# To avoid index calculations, maintain i = INDEX(c), nc = this.corners[c].opp, ni = INDEX(nc)
	while (DEAD(this, SIMPLEX(nc)))
		ni::Int32 = INDEX(nc)
		@inbounds off::Int32 = indoff[ni+1, j+1, i+1] 	# where j goes relative to i is our new i
		j = ni
		i = ni + off 	# fix new j, i, c, and try neighbor
		c = nc + off
		@inbounds nc = this.corners.opp[c]
	end
		
	@inbounds nc = this.corners.opp[nc] 	# go to new simplex
	@inbounds @assert this.corners.vertex[nc] == pv "Expected to find new simplex using pv after walking dead simplices"

	@inbounds setCornerVC!(this, newb, v, Int32(nc - 3 + invdrop[i+1, j+1]-1))

	return v
end

# ╔═╡ e27b7858-ec7e-493d-9979-87f1b96c792e
md"### Insert!"

# ╔═╡ c53bd236-331e-47af-8be8-59348076d145
# Testing insert
# begin
# 	this::StateType = StateType()
# 	nVert::Int32 = 7
# 	vertArray::Vector{Point} = Vector{Point}(undef, nVert)
# 	# vertArray[1] = Point(0, 0, 0, 1)
# 	# vertArray[2] = Point(3, 2, 2)
# 	# vertArray[3] = Point(1, 2, 1)
# 	# vertArray[4] = Point(2, 3, 0)
# 	# vertArray[5] = Point(3, 1, 1)
# 	# vertArray[6] = Point(2, 2, 1)

# 	vertArray[1] = Point(0, 0, 0, 1)
# 	vertArray[2] = Point(0, 0, 0)
# 	vertArray[3] = Point(2, 0, 0)
# 	vertArray[4] = Point(0, 2, 0)
# 	vertArray[5] = Point(0, 0, 2)
# 	vertArray[6] = Point(1, -2, 1)
# 	vertArray[7] = Point(5, 5, 5)

# 	# vertArray[1] = Point(0, 0, 0, 1)
# 	# vertArray[2] = Point(-3, 0, 0)
# 	# vertArray[3] = Point(3, 0, 0)
# 	# vertArray[4] = Point(0, 3, 0)
# 	# vertArray[5] = Point(0, 0, 3)
# 	# vertArray[6] = Point(0, 1, 1)
# 	# vertArray[7] = Point(5, 5, 5)

# 	d3Initialize!(this, vertArray, nVert)
# 	# this
	
# 	index::Int32 = 6
# 	p = d3LocSphere(this, this.verts[index], this.liveSimp)[1]
# 	d3Insert!(this, index, p)

# 	index = 7
# 	p = d3LocSphere(this, this.verts[index], this.liveSimp)[1]
# 	d3Insert!(this, index, p)

# 	auditCorners(this, Int32(1))
# end

# ╔═╡ 59012fdf-efa9-4537-914f-907866c43d31
md"# Compact Corners"

# ╔═╡ 78dc7bf0-105b-444d-889e-3cf2be379840
# TODO?

# ╔═╡ fa8e47a1-0f7e-4f62-be54-e15002eaefdc
md"# Audit"

# ╔═╡ 7063fb91-05f8-456d-b167-ec6122e2a6d1
# Print corners
begin
	function cornerPrint(this::StateType, c::Int32)
		println("Corner: ", c, ", Simplex: ", SIMPLEX(c), ", Index: ", INDEX(c))

		@inbounds println("( ", 
			infiniteV(this.corners.vertex[c], this.corners.vertex[4]) ? 0 : 1, " ",
			this.corners.vertex[c].x, " ",
			this.corners.vertex[c].y, " ",
			this.corners.vertex[c].z, " ",
			this.corners.vertex[c].sq, " ) opp:",
			this.corners.opp[c], "(",
			SIMPLEX(this.corners.opp[c]), ",", 
			INDEX(this.corners.opp[c]), ")"
		)
	end;

	function cornerPrint4(this::StateType, c::Int32)
		b::Int32 = BASECORNER(c)
		if (@inbounds this.sphs[1] != undef) 
			sp::Sphere = this.sphs[SIMPLEX(b)]
			println("disp('Sphere(", SIMPLEX(b), ") = < ",
				sp.x, " ",
				sp.y, " ",
				sp.z, " ",
				sp.sq, " >')")
		end

		println("DetCheckH([")
		for k in 0:3
			@inbounds println(" ", 
				infiniteV(this.corners.vertex[b+k], this.corners.vertex[4]) ? 0 : 1, " ",
				this.corners.vertex[b+k].x, " ",
				this.corners.vertex[b+k].y, " ",
				this.corners.vertex[b+k].z, " ",
				this.corners.vertex[b+k].sq, ";")
		end
		println("]);")

		for j::Int32 in 0:3
			cornerPrint(this, b+j)
		end
	end;
end

# ╔═╡ b4142133-c7d6-4ec8-a709-d69cdc90be82
#=
Inserts the point this.vert[vi] that is contained in sphere p into the delaunay triangulation stored in this. 
(d3LocSphere may be used to optain p.)
=#

function d3Insert!(this::StateType, vi::Int32, p::Int32)
	pv::Point = this.verts[vi]
	
	#=
	Simplices containing pv are "dead", and are pushed onto kill stack.
	We use DFS with stack to find them and kill them
	At live-dead boundary, we save dead simplices on stack nhbr,
		then make new simplices and hook in to live by setting the last opp pointer.

	Invariants/operations: 
		Simplex p is marked alive or dead on first visit.
		Corner c is pushed on stack when SIMPLEX(this.corners.opp[c]) is marked dead

	On termination, stack nhbr contains dead corners with live neighbors
		that have new simplices (so this.corners.opp[nhbr] != this.corners.opp[this.corners.opp[nhbr]] temporarily.)
	Stack kill contains old simplices for final recycling
	=#

	# Initialize stacks
	this.dfs = Stack() 	# DFS stack holds corners opposite dead simplices
	this.idfs = Stack() 	# iDFS stack holds corners opposite infinite simplices with pv on boundary
						# 	(these are special case: dead, but don't propagate)
	this.nhbr = Stack() 		# stack for dead corners with live nhbr simplices
	this.kill = Stack() 		# stack of dead simplices to recycle

	b::Int32 = CORNER(p, Int32(0))
	push!(this.kill, p) 	# Kill simp initial p
	KILL(this, p)
	@inbounds push!(this.dfs, this.corners.opp[b]) 	# stack neighbors
	b+=1
	@inbounds push!(this.dfs, this.corners.opp[b])
	b+=1
	@inbounds push!(this.dfs, this.corners.opp[b])
	b+=1
	@inbounds push!(this.dfs, this.corners.opp[b])

	# Perform the DFS
	DFS!(this, pv, p)

	# Check the neighbors of infinite simps
	while (!isEmpty(this.idfs))
		c::Int32 = pop!(this.idfs)
		p = SIMPLEX(c)

		if (!DEAD(this, SIMPLEX(this.corners.opp[c])))
			println()
			println("Problem!:")
			for i::Int32 in 1:this.maxSimps
				println(vi)
				if (this.active[i] > 0) 
					cornerPrint4(this, Int32(i*4))
				end
			end	
		end
		
		@inbounds @assert DEAD(this, SIMPLEX(this.corners.opp[c])) "dfs stack element w/ non-dead neighbor"

		if (DEAD(this, p)) 	# dead already
			continue
		end

		@inbounds @assert DEAD(this, SIMPLEX(this.corners.opp[c])) "Live corner c should have dead neighbor"

		@inbounds push!(this.nhbr, this.corners.opp[c])	# remember old corner to hook simplex into mesh later
		STARTSIMP(this) 	# make new simplex liveSimp
		newb = CORNER(this.liveSimp, Int32(3)) 	# Last corner of new tetra
		setCornerVCN!(this, newb, pv, c) 	# last corner is pv; also set oppositie corner c. Do rest later.
	end

	# Now we have stack of dead neighbors of live simps, and we've hooked new simps to them
	while (!isEmpty(this.nhbr))
		dead::Int32 = pop!(this.nhbr) 	# dead simplex
		jdead::Int32 = INDEX(dead) 		# index of dropped corner

		@assert DEAD(this, SIMPLEX(dead)) "corner on nhbr stack is not dead!?"
		
		@inbounds newb::Int32 = this.corners.opp[this.corners.opp[dead]] - 3 	# base of new tetra
												# This was set above

		dead -= jdead 	# just use base of dead one

		# new tetra has 0, 1, 2, 3 = pv 
		# corresponding old indices before jdead is dropped:
		# 	drop[j][0], ..., drop[j][3], (no corresp to pv)

		# Set the corners of the new simplex

		# v0::Point = newSimpCorner!(this, pv, newb, dead, 1)
		# newb += 1
		# v1::Point = newSimpCorner!(this, pv, newb, dead, 2)
		# newb += 1
		# v2::Point = newSimpCorner!(this, pv, newb, dead, 3)
		# newb += 1
		
		
		### v0 ###
		# TODO: Maybe fix the drop table values to avoid the unnecessary -1
		# TODO: Also fix indices that have +1's 
		j::Int32 = jdead
		@inbounds i::Int32 = drop[j+1, 1]-1 	# old index of new corner 1
		c::Int32 = dead + i 		# Note: i = INDEX(c)
		@inbounds v0::Point = this.corners.vertex[c] 	# copy vertex v0
		@inbounds nc::Int32 = this.corners.opp[c] 		# go to neighbor

		# In simplex opp c, find new location of j. That is new c. New j = INDEX(c.opp) + 1.
		# To avoid index calculations, maintain i = INDEX(c), nc = this.corners[c].opp, ni = INDEX(nc)
		while (DEAD(this, SIMPLEX(nc)))
			ni::Int32 = INDEX(nc)
			@inbounds off::Int32 = indoff[ni+1, j+1, i+1] 	# where j goes relative to i is our new i
			j = ni
			i = ni + off 	# fix new j, i, c, and try neighbor
			c = nc + off
			@inbounds nc = this.corners.opp[c]
		end
			
		@inbounds nc = this.corners.opp[nc] 	# go to new simplex
		@inbounds @assert this.corners.vertex[nc] == pv "Expected to find new simplex using pv after walking dead simplices"

		@inbounds setCornerVC!(this, newb, v0, Int32(nc - 3 + invdrop[i+1, j+1]-1))
		newb += 1

		

		### v1 ###
		j = jdead
		@inbounds i = drop[j+1, 2]-1 	# old index of new corner 1
		c = dead + i 		# Note: i = INDEX(c)
		@inbounds v1::Point = this.corners.vertex[c] 	# copy vertex v1
		@inbounds nc = this.corners.opp[c] 		# go to neighbor

		# In simplex opp c, find new location of j. That is new c. New j = INDEX(c.opp) + 1.
		# To avoid index calculations, maintain i = INDEX(c), nc = this.corners.opp[c], ni = INDEX(nc)
		while (DEAD(this, SIMPLEX(nc)))
			ni::Int32 = INDEX(nc)
			@inbounds off::Int32 = indoff[ni+1, j+1, i+1] 	# where j goes relative to i is our new i
			j = ni
			i = ni + off 	# fix new j, i, c, and try neighbor
			c = nc + off
			@inbounds nc = this.corners.opp[c]
		end
			
		@inbounds nc = this.corners.opp[nc] 	# go to new simplex
		@inbounds @assert this.corners.vertex[nc] == pv "Expected to find new simplex using pv after walking dead simplices"

		@inbounds setCornerVC!(this, newb, v1, Int32(nc - 3 + invdrop[i+1, j+1]-1))
		newb += 1

		

		### v2 ###
		j = jdead
		@inbounds i = drop[j+1, 3]-1 	# old index of new corner 1
		c = dead + i 		# Note: i = INDEX(c)
		@inbounds v2::Point = this.corners.vertex[c] 	# copy vertex v2
		@inbounds nc = this.corners.opp[c] 		# go to neighbor

		# In simplex opp c, find new location of j. That is new c. New j = INDEX(c.opp) + 1.
		# To avoid index calculations, maintain i = INDEX(c), nc = this.corners.opp[c], ni = INDEX(nc)
		while (DEAD(this, SIMPLEX(nc)))
			ni::Int32 = INDEX(nc)
			@inbounds off::Int32 = indoff[ni+1, j+1, i+1] 	# where j goes relative to i is our new i
			j = ni
			i = ni + off 	# fix new j, i, c, and try neighbor
			c = nc + off
			@inbounds nc = this.corners.opp[c]
		end
			
		@inbounds nc = this.corners.opp[nc] 	# go to new simplex
		@inbounds @assert this.corners.vertex[nc] == pv "Expected to find new simplex using pv after walking dead simplices"

		@inbounds setCornerVC!(this, newb, v2, Int32(nc - 3 + invdrop[i+1, j+1]-1))
		newb += 1 # I think this line is unnecessary?


		@inbounds c = this.corners.opp[this.corners.opp[dead+jdead]]

		if (v0 != this.corners.vertex[CORNERINOPP(this, 0, c)])
			println(v0)
			println(this.corners.vertex[CORNERINOPP(this, 0, c)])
		end

		@inbounds @assert v0==this.corners.vertex[CORNERINOPP(this, 0, c)] "v0 does not line up"
		@inbounds @assert v1==this.corners.vertex[CORNERINOPP(this, 1, c)] "v1 does not line up"
		@inbounds @assert v2==this.corners.vertex[CORNERINOPP(this, 2, c)] "v2 does not line up"

		@inbounds this.sphs[SIMPLEX(newb)] = makeSphere(v0, v1, v2, pv, this.verts[1])
	end
	
end

# ╔═╡ a22b166b-6c53-4858-a99d-fbfb37744443
# Audit corners
function auditCorners(this::StateType, sphereCheck::Int32)
	for p::Int32 in 1:this.maxSimps
		if (DEAD(this, p)) 
			# Don't audit dead simps
			continue
		end

		b::Int32 = CORNER(p, Int32(0))

		# Per corner checks
		for c in b:CORNER(p, Int32(3))
			# Check opposite
			@inbounds i::Int32 = this.corners.opp[c]
			if (@inbounds this.corners.opp[i] != c) 
				println("AUDIT: wrong opp.opp")
				cornerPrint(this, c)
				cornerPrint(this, i)
			end

			# Make sure it doesn't point to the same vertex
			if (@inbounds this.corners.vertex[c] == this.corners.vertex[i])
				println("AUDIT: Same vertex  ", 
					"corners.vertex[", c, "(", SIMPLEX(c), ",", INDEX(c), ")] == ",
					"opp.vertex[", i, "(", SIMPLEX(i), ",", INDEX(i), ")]"
					)
				cornerPrint4(this, c)
				cornerPrint4(this, i)
			end

			# Check CORNERINOPP table
			for j in 0:3
				if (j != INDEX(c) && CORNERINOPP(this, j, c) < 3)
					if (@inbounds this.corners.vertex[BASECORNER(c)+j] != this.corners.vertex[CORNERINOPP(this, j, c)])
						@inbounds println("AUDIT: Bad auditCornerInOpp(", j, ",", c, ") = ", CORNERINOPP(this, j, c),
							"  since vertex ",
							this.corners.vertex[BASECORNER(c)+j], " != ",
							this.corners.vertex[CORNERINOPP(j, c)]
							)
						cornerPrint4(this, BASECORNER(c))
						cornerPrint4(this, BASECORNER(CORNERINOPP(this, j, c)))
						break
					end
				end
			end

			for j in 0:3
				k::Int32 = CORNERINOPP(this, j, c)
				if (@inbounds SIMPLEX(k) != SIMPLEX(this.corners.opp[c]) || this.corners.vertex[b+j] != this.corners.vertex[k])
					if (@inbounds SIMPLEX(k) != SIMPLEX(this.corners.opp[c]))
						@inbounds println("AUDIT: CORNERINOPP(", j, ",", c, ") ==> ",
							k - this.corners.opp[c], ": Accessing [",
							c, ",", INDEX(c), "][",
							this.corners.opp[c], INDEX(this.corners.opp[c]), "][",j, "]"
						)
					else
						if (b + j == c)
							if (@inbounds k == this.corners.opp[c])
								continue 	# These vertices are supposed to differ; don't flag them
							else
								println("AUDIT: CORNERINOPP(", j, ",", c, ") says ",
									b+j, "(", SIMPLEX(b), ",", j, ") and ",
									k,   "(", SIMPLEX(k), ",", INDEX(k), ") shouldn't happen"
								)
							end
						else 
							println("b: ", b)
							println("j: ", j)
							println("c: ", c)
							println("AUDIT: CORNERINOPP(", j, ",", c, ") says ",
								b+j, "(", SIMPLEX(b), ",", j, ") and ",
								k,   "(", SIMPLEX(k), ",", INDEX(k), ") should agree"
								)
							cornerPrint4(this, c)
							cornerPrint4(this, Int32(b+j))
							cornerPrint4(this, k)

							# The 'j' in this assert was an 'i' in the C code, but that did not seem correct, as this value should be between 1 and 4
							@inbounds @assert SIMPLEX(this.corners.opp[c]) == SIMPLEX(CORNERINOPP(this, Int64(j), c)) "CORNERINOPP screws up simplices"
						end
					end

					cornerPrint4(this, b)
					cornerPrint4(this, BASECORNER(k))
					break
				end
			end

			# Check Sphere opposite corner c
			if (sphereCheck > 0) 	# Check sphere opposite corner
				@inbounds k = this.corners.opp[c]
				@inbounds sp = this.sphs[SIMPLEX(k)]
				d = 0
				if (@inbounds infiniteV(this.corners.vertex[c], this.verts[1]))
					d = sp.sq
				else
					@inbounds vp = this.corners.vertex[CORNER(SIMPLEX(k), Int32(3))]
					@inbounds d = spdot(sp, this.corners.vertex[c], vp)
				end

				if (d < 0)
					@inbounds println("disp('AUDIT: corner ", c, " vertex ( ",
						this.corners.vertex[c].x, " ",
						this.corners.vertex[c].y, " ",
						this.corners.vertex[c].z, " ",
						this.corners.vertex[c].sq, " ) in sphere ", 
						SIMPLEX(k), "(", k, ") = ", d, "');"
					)
					cornerPrint(this, c)
					cornerPrint(this, k)
				end
			end
		end

		if (sphereCheck > 0) 	# Check sphere sqs (orient dets)
			@inbounds sp = this.sphs[p] 	# only sphers using point @ infinity have sq==0; none have sq<0
			b = CORNER(p, Int32(0))

			if (@inbounds sp.sq < 0 || (sp.sq == 0 && !infiniteV(this.corners.vertex[b], this.verts[1]) && !infiniteV(this.corners.vertex[b+1], this.verts[1])))
				println("disp('AUDIT: sq<=0 in simplex ", CORNER(p,Int32(0)), "(", p, ") = ",
					sp.sq, "');"					
				)
				cornerPrint4(this, b)
			end
		end
	end
end

# ╔═╡ 4f266903-8c4f-4c57-b5a6-8f629b05c0c6
md"# Batch"

# ╔═╡ 87b13c4d-416e-4a78-a182-82a0bd9130f7
# Utility functions for generating data
begin
	# Generate a random value from 0 to one under the highest coordinate allowed above
	@inline function randomCoord()::Int32
		return rand(Int32) & COORDMASK
	end
	
	# Generate a random Point
	@inline function randomPoint()::Point
		x::Int32 = randomCoord()
		y::Int32 = randomCoord()
		z::Int32 = randomCoord()
		return Point(x, y, z)
	end
	
	# Genereate 'num' random points
	# Adds point @ infinity in the first index
	function generateRandomPoints(num::Int32)::Vector{Point}
		output::Vector{Point} = Vector{Point}(undef, num)

		for i in 1:num
			@inbounds output[i] = randomPoint()
		end

		output[1] = Point(0, 0, 0, 1) 	# Point @ infinity must come first
		return output
	end

	# Generate points in a num by num by num grid
	# Adds point @ infinity in the first index
	function generateGridPoints(num::Int32)::Vector{Point}
		output::Vector{Point} = Vector{Point}(undef, num*num*num+1)

		for x in 0:num-1
			for y in 0:num-1
				for z in 0:num-1
					@inbounds output[z+y*num+x*num*num+2] = Point(x, y, z)
				end
			end
		end

		output[1] = Point(0, 0, 0, 1) 	# Point @ infinity must come first
		return output
	end

	# Generate points (t, t^2, t^3) to create skinny simplices
	function generateCurvePoints(num::Int32)::Vector{Point}
		output::Vector{Point} = Vector{Point}(undef, num+1)

		for t in 0:num-1
			output[t+2] = Point(t, t*t, t*t*t)
		end

		output[1] = Point(0, 0, 0, 1) 	# Point @ infinity must come first
		return output
	end
end;

# ╔═╡ 6420ebbb-e126-45d4-b894-30075c73ef4a
#=
d3batch takes lifted input vertices and returns a (compact) corner table for Delaunay.
REQUIRES that the first point is @ infinity,
and that the first 5 points are in general position. (I should verify or relax this.)
Since all points have radii assigned already, it can compute power diagrams.
=#

# Input vertices and number of vertices
# Output corner table and number of corners
function d3Batch(vertArray::Vector{Point}, nVert::Int32)
	itry::Int32 = 5
	j::Int32 = 0
	
	this::StateType = StateType();
	while (!d3Initialize!(this, vertArray, nVert)) 
		println("Permute points to get the first five vertices to be in general position")

		j = (j+1)%4 + 2
		vertArray[j], vertArray[itry] = vertArray[itry], vertArray[j]
		itry+=1

		if (itry == nVert) 
			println("Can't find vertices to initialize. Quit")
			return ([], -1) 	# Empty return to symbolize failure
		end
	end

	for i in 1:nVert
		println("$(i) $(vertArray[i])")
	end

	maxLocate = locateSideCnt = sphereCnt = startTetraCnt = freeTetraCnt = inSphereCnt = 0

	for vi::Int32 in 5:nVert
		# incrementally insert vert[vi]
		# LOCATE: find some tetrahedron with sphere strictly containing vert[vi]
		(p,k) = d3LocSphere(this, this.verts[vi], this.liveSimp)
		if ( k < 1 )
			continue
		end

		# Insert vertex vi, which is in sphere of simplex p
		for i::Int32 in 1:this.maxSimps
			println("Before ", vi)
			if (this.active[i] > 0) 
				cornerPrint4(this, Int32(i*4))
			end
		end	
		
		d3Insert!(this, vi, p)

		while (!isEmpty(this.kill))
			# Recycle memory of dead simplices
			p = pop!(this.kill)
			FREESIMP(this, p)
		end
		
			
	end

	# Count the corners
	nSimps = 0
	for j::Int32 in 1:this.maxSimps
		this.active[j] = DEAD(this, j) ? -1 : (nSimps += 1)
	end

	# Add compact corners?

	# Audit the corners
	auditCorners(this, Int32(1))

	# Return number of corners and 'this'
	return (nSimps, this)

	# Return just the corners
	return (this.corners, nSimps)
end

# ╔═╡ 3813e5e7-1995-492f-98a7-01d424f25d82
md"# Testing"

# ╔═╡ bb33f4e7-2587-4e0e-93af-0ac4f580bba8
# begin
# 	### batch testing ###
# 	# Manually entered points
	
# 	nVert::Int32 = 7
	
# 	vertArray::Vector{Point} = Vector{Point}(undef, nVert)

# 	vertArray[1] = Point(0, 0, 0, 1) 	# Point @ infinity must come first

# 	# Start with some known points to test.
# 	vertArray[2] = Point(0, 0, 0)
# 	vertArray[3] = Point(2, 0, 0)
# 	vertArray[4] = Point(0, 2, 0)
# 	vertArray[5] = Point(0, 0, 2)
# 	vertArray[6] = Point(1, -2, 1) 
# 	vertArray[7] = Point(2, 2, 2)

# 	d3Batch(vertArray, nVert)
# end

# ╔═╡ 62fbbd77-5f7b-4027-8f3f-67ed80765056
begin
	### batch testing ###
	# Manually entered points
	
	nVert::Int32 = 13
	
	vertArray::Vector{Point} = Vector{Point}(undef, nVert)

	vertArray[1] = Point(0, 0, 0, 1) 	# Point @ infinity must come first

	# Start with some known points to test.
	vertArray[2] = Point(0, 0, 0)
	vertArray[3] = Point(0, 0, 1)
	vertArray[4] = Point(0, 1, 0)
	vertArray[5] = Point(1, 0, 0)
	
	vertArray[6] = Point(0, 0, 2)
	
	vertArray[7] = Point(0, 1, 1) 
	vertArray[8] = Point(0, 1, 2)

	vertArray[9]  = Point(0, 2, 1) 
	vertArray[10] = Point(0, 2, 0)
	vertArray[11] = Point(0, 2, 2)
	
	
	vertArray[12] = Point(1, 0, 1) 
	vertArray[13] = Point(1, 0, 2)

	(nSimps, this) = d3Batch(vertArray, nVert)
	
end

# ╔═╡ 74febf8f-ecdc-455c-9d04-6024b89f2e6a
begin
	### batch testing ###
	# Random points

	# random_nVert::Int32 = 3000
	
	# random_vertArray = generateRandomPoints(random_nVert)

	# d3Batch(random_vertArray, random_nVert)
end

# ╔═╡ 9f72cd6f-68ed-4661-bf7c-efa0ecfa99fb
begin
	### batch testing ###
	# Grid points

	# gridSize::Int32 = 3
	
	# grid_vertArray = generateGridPoints(gridSize)
	
	# d3Batch(grid_vertArray, Int32(13)) #gridSize*gridSize*gridSize)
end

# ╔═╡ 35174515-ba51-4230-8117-84a4884fbb1f


# ╔═╡ 1b00e7db-1e9c-4969-b078-8b75adad721c
begin
	### batch testing ###
	# Curve points (t, t^2, t^3)

	# curve_nVert::Int32 = 100
	
	# curve_vertArray = generateCurvePoints(curve_nVert)
	
	# d3Batch(curve_vertArray, curve_nVert)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[compat]
BenchmarkTools = "~1.3.1"
StructArrays = "~0.6.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.1"
manifest_format = "2.0"
project_hash = "f2209c5b85cf8a0221b6b1295d1cad8c0d0cd709"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.DataAPI]]
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "13237798b407150a6d2e2bce5d793d7d9576e99e"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.13"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
# ╟─3138335c-5ce5-4f70-bf3e-7d093ce8de13
# ╠═93e62420-39d1-4813-9bd5-7fa12286d190
# ╟─9f95fab7-f3d3-46c1-ac32-537b84196a35
# ╠═23f069f2-c3a8-47ff-8af5-1f3a242e6bb4
# ╟─512b21a0-1e6f-44b8-a05c-38d45ec7467f
# ╠═cada11d4-0aa9-4b74-a2ff-593b2497e9f6
# ╠═083f724b-2009-4a49-ae1e-61363c3d7165
# ╠═2e9a743d-78d2-4a89-a400-3f6ade10b7e0
# ╟─c368ce56-b068-45ac-979f-084e3175c8e0
# ╟─75695a4f-7c0c-4879-88d8-c2fc7acbccca
# ╠═52f11677-636c-4a90-9eeb-cf5e4b948bfe
# ╠═81ce6604-91ef-4751-b4ff-1535dad19c40
# ╠═77719559-7da6-49ac-b5f3-e5de2428bb8a
# ╟─f92bf3ea-b573-4790-8b84-348ddef64aca
# ╠═0344f8c1-49d9-482b-b9e0-60ca3bfb8bab
# ╠═8ddc46f4-2e4b-4184-9a49-5c5d6d3be112
# ╟─bf22c408-6d8d-45d3-a4d7-cbda3dd9fa1a
# ╟─8ced7d59-909a-46b1-954e-0c91324ec2d1
# ╟─9db9a542-ed11-4408-bbe9-a2b1d56bca96
# ╟─bd61ed12-f3c8-452d-b848-b4d9994164af
# ╠═c21855ca-c54d-4e4b-b68c-b3b548eda1db
# ╠═d52396af-0dfb-4fa0-9942-fba52d282eae
# ╟─27daafcc-f7b1-43bc-8689-70fef984fa5f
# ╠═3f8be73a-4ce9-4ef1-bb94-e572e82ac382
# ╟─08e8aeb5-8d8e-44b7-bdd3-1f8e036611b6
# ╠═f4fff776-1f97-4bd9-9fd5-7d22780a7028
# ╠═25fb3749-b637-4095-b277-8b3325161a25
# ╟─d7cdc261-232c-4560-95f1-52aa2e4074c7
# ╠═894587ef-8eee-4eeb-8354-2dbbda06a852
# ╟─492d064a-fc8e-45d6-a548-e206d7a6e65b
# ╠═b112c2d3-a38e-41e0-a2be-839faa770b8e
# ╟─9c700c27-565f-4a9b-b031-229f434d0512
# ╟─5313f41e-19d7-4257-82e0-d04db110050b
# ╟─f5b545bc-ab45-4c6c-9843-c076b79d4915
# ╟─f3d6901e-8777-4505-91c6-4b9be755f5c3
# ╟─9f387d99-7848-4bbb-83e3-160f6fe08dfb
# ╟─f9b13d44-ecde-4133-86de-bdb9506dd99a
# ╠═4de76870-e6c6-48de-a954-75e751ef2944
# ╟─fcf8c6cd-d50d-4991-a5dd-5fc0bcc8483b
# ╠═dd5593da-8154-4988-820f-a73a24c9a3d6
# ╟─e27b7858-ec7e-493d-9979-87f1b96c792e
# ╠═b4142133-c7d6-4ec8-a709-d69cdc90be82
# ╟─c53bd236-331e-47af-8be8-59348076d145
# ╟─59012fdf-efa9-4537-914f-907866c43d31
# ╠═78dc7bf0-105b-444d-889e-3cf2be379840
# ╟─fa8e47a1-0f7e-4f62-be54-e15002eaefdc
# ╠═7063fb91-05f8-456d-b167-ec6122e2a6d1
# ╟─a22b166b-6c53-4858-a99d-fbfb37744443
# ╟─4f266903-8c4f-4c57-b5a6-8f629b05c0c6
# ╠═87b13c4d-416e-4a78-a182-82a0bd9130f7
# ╠═6420ebbb-e126-45d4-b894-30075c73ef4a
# ╟─3813e5e7-1995-492f-98a7-01d424f25d82
# ╠═bb33f4e7-2587-4e0e-93af-0ac4f580bba8
# ╠═62fbbd77-5f7b-4027-8f3f-67ed80765056
# ╠═74febf8f-ecdc-455c-9d04-6024b89f2e6a
# ╠═9f72cd6f-68ed-4661-bf7c-efa0ecfa99fb
# ╠═35174515-ba51-4230-8117-84a4884fbb1f
# ╠═1b00e7db-1e9c-4969-b078-8b75adad721c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
