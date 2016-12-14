#########################################################################################
#
# napl_layer.jl - a julia-language script to model the movement of a NAPL layer floating 
# on top if a single-layer, unconfined aquifer
#
# volume elements are defined for a single 2-D layer (can be unstructured)
# a hidden, overlying second layer models the NAPL
#
#########################################################################################

type Monitor 					# monitoring well "data logger"
	time::Float64
	hw::Float64
	hn::Float64
end


type Fluid
	name::AbstractString
	rho::Float64
end


type Node
	x::Float64
	y::Float64
	z0::Float64
	area::Float64 									# planar area of volume element (fluid volume = this number * fluid thickness * Sy)
	Kw::Float64 									# effective hydraulic conductivities (water and NAPL)
	Kn::Float64
	Sy::Float64 									# specific storage
	hw::Float64 									# water and NAPL thickness
	hn::Float64
	Qw::Float64 									# volumetric flux source term for water and NAPL (constant in time)
	Qn::Float64 
	connect_list::Array{Tuple{Int64, Int64}} 		# list of connected elements/nodes (connection number, connecting node index number)
	sigma_w::Float64 								# water and NAPL total conductance (updated each time step)
	sigma_n::Float64
end


type Connect
	node_1::Int64 				# index numbers for connected nodes
	node_2::Int64
	dx::Float64					# inter-node distance (for computing gradient)
	len_inf::Float64 			# volume element connection interfacial length
	conduct_w::Float64 			# connection conductance (water and NAPL); K*A/dx
	conduct_n::Float64
end


type Knobs
	gamma::Float64 				# time-stepping weighting factor for implicit solution scheme
	dt_init::Float64 			# initial, minimum, and maximum time step size
	dt_min::Float64
	dt_max::Float64
	dh_max::Float64 			# maximum change in fluid thickness (water or NAPL), per time step
	dt_decrease::Float64		# time step reduction and increase factors
	dt_increase::Float64
end


function GetKnobs()
	# read numerical model "knobs" from file
	data = readdlm("knobs.txt", '\t', header=false)
	gamma = Float64(data[1, 2])
	dt_init = Float64(data[2, 2])
	dt_min = Float64(data[3, 2])
	dt_max = Float64(data[4, 2])
	dh_max = Float64(data[5, 2])
	dt_decrease = Float64(data[6, 2])
	dt_increase = Float64(data[7, 2])
	knobs = Knobs(gamma, dt_init, dt_min, dt_max, dh_max, dt_decrease, dt_increase)
	println("Read in computational knobs.")
	return knobs
end


function ReadFluids()
	# read in fluid properties (density)
	fluid = Fluid[]
	data = readdlm("fluids.txt", '\t', header=true)
	for i = 1:2
		name = data[1][i, 1]
		rho = Float64(data[1][i, 2])
		push!(fluid, Fluid(name, rho))
	end
	println("Read fluid properties.")
	return fluid
end


function ReadNodes()

	# read in nodes file and populate node type array
	node = Node[]
	num_nodes = 0
	data = readdlm("nodes.txt", '\t', header=true)
	for i = 1:size(data[1], 1)
		x = Float64(data[1][i, 1])
		y = Float64(data[1][i, 2])
		z0 = Float64(data[1][i, 3])
		num = Int64(data[1][i, 4])
		x_step = Float64(data[1][i, 5])	
		y_step = Float64(data[1][i, 6])	
		z_step = Float64(data[1][i, 7])		
		area = Float64(data[1][i, 8])	
		hw = Float64(data[1][i, 9])	
		hn = Float64(data[1][i, 10])			
		Kw = Float64(data[1][i, 11])	
		Kn = Float64(data[1][i, 12])			
		Sy = Float64(data[1][i, 13])					
		Qw = Float64(data[1][i, 14])	
		Qn = Float64(data[1][i, 15])				
		for j = 1:num
			# for each like node; connecting node list and sigma terms will be created later from connections		
			push!(node, Node(x+j*x_step, y+j*y_step, z0+j*z_step, area, Kw, Kn, Sy, hw, hn, Qw, Qn, Tuple{Int64, Int64}[], 0., 0.))
			num_nodes += 1
		end
	end
	
	# write to node summary file ...
	fname = "nodes_summary.csv"
	csvfile = open(fname,"w")
	line_out = "node" * "," * "x" * "," * "y" * "," * "z0" * "," * "area" * "," * "Kw" * "," * "Kn" * "," * "Sy" *
		"," * "hw" * "," * "hn" * "," * "Qw" * "," * "Qn"
	println(csvfile, line_out)	
	for i = 1:num_nodes
		line_out = string(i) * "," * string(node[i].x) * "," * string(node[i].y) * "," * string(node[i].z0) * "," * string(node[i].area) *
			"," * string(node[i].Kw) * "," * string(node[i].Kn) * "," * string(node[i].Sy) * "," * string(node[i].hw) * "," * string(node[i].hn) *
			"," * string(node[i].Qw) * "," * string(node[i].Qn)
		println(csvfile,line_out)
	end
	close(csvfile)		
	println("Processed nodes.")	
	
	return node, num_nodes

end


function Conduct(node::Array{Node,1}, i::Int64, j::Int64, dx::Float64, len_inf::Float64)
	# compute/update conductance between two nodes (per fluid)
	conduct_w = mean([node[i].Kw * node[i].hw, node[j].Kw * node[j].hw]) * len_inf / dx
	conduct_n = mean([node[i].Kn * node[i].hn, node[j].Kn * node[j].hn]) * len_inf / dx	
	return conduct_w, conduct_n
end


function ReadConnects(node::Array{Node,1})

	# read in connections file and populate node type array
	connect = Connect[]
	num_connects = 0
	data = readdlm("connects.txt", '\t', header=true)
	for i = 1:size(data[1], 1)
		node_1 = Int64(data[1][i, 1])
		node_2 = Int64(data[1][i, 2])
		num = Int64(data[1][i, 3])
		step_1 = Int64(data[1][i, 4])	
		step_2 = Int64(data[1][i, 5])	
		delta = Float64(data[1][i, 6])		
		len_inf = Float64(data[1][i, 7])	
		for j = 1:num
			conduct_w, conduct_n = Conduct(node, node_1, node_2, delta, len_inf)
			push!(connect, Connect(node_1+j*step_1, node_2+j*step_2, delta, len_inf, conduct_w, conduct_n)) 	# for each like connection
			num_connects += 1
		end
	end
	
	# write to connection summary file ...
	fname = "connects_summary.csv"
	csvfile = open(fname,"w")
	line_out = "connection" * "," * "node_1" * "," * "node_2" * "," * "delta" * "," * "len_inf"
	println(csvfile, line_out)	
	for i = 1:num_connects
		line_out = string(i) * "," * string(connect[i].node_1) * "," * string(connect[i].node_2) *
			"," * string(connect[i].dx) * "," * string(connect[i].len_inf)
		println(csvfile,line_out)
	end
	close(csvfile)		
	println("Processed connections.")	
	
	return connect, num_connects

end

function LHSmatrix(connect::Array{Connect,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64, num_nodes::Int64, alpha::Float64)

	# fill out the LHS of the equation matrix and record row-column index positions

	row_index = Int64[] 					# indexing system for sparse matrix
	col_index = Int64[]
	data = Float64[]
	
	# diagonal elements
	for (i, nd) in enumerate(node) 
		push!(row_index, i) 								
		push!(col_index, i)				
		push!(data, nd.area*nd.Sy/dt + knobs.gamma*nd.sigma_w)		# water-water
		push!(row_index, i) 								
		push!(col_index, i + num_nodes)				
		push!(data, knobs.gamma * alpha * nd.sigma_w) 				# water-NAPL
		push!(row_index, i + num_nodes) 								
		push!(col_index, i)				
		push!(data, knobs.gamma * nd.sigma_n)	 					# NAPL-water	
		push!(row_index, i + num_nodes) 								
		push!(col_index, i + num_nodes)				
		push!(data, nd.area*nd.Sy/dt + knobs.gamma*nd.sigma_n) 		# NAPL-NAPL
	end	
	
	# non-diagonal elements
	for (i, nd) in enumerate(node) 
		for cn in nd.connect_list
			push!(row_index, i) 								
			push!(col_index, cn[2])
			push!(data, -knobs.gamma * connect[cn[1]].conduct_w) 			# water-water
			push!(row_index, i) 								
			push!(col_index, cn[2] + num_nodes)
			push!(data, -knobs.gamma * alpha * connect[cn[1]].conduct_w) 	# water-NAPL	
			push!(row_index, i + num_nodes) 
			push!(col_index, cn[2])				
			push!(data, -knobs.gamma * connect[cn[1]].conduct_n) 			# NAPL-water	
			push!(row_index, i + num_nodes) 
			push!(col_index, cn[2] + num_nodes)	
			push!(data, -knobs.gamma * connect[cn[1]].conduct_n) 			# NAPL-NAPL
			
		end
	end
	
	return data, row_index, col_index
	
end


function UpdateLHS(connect::Array{Connect,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64, num_nodes::Int64, alpha::Float64)

	# fill out the LHS of the equation matrix (row-column index positions already have been recorded)

	data = Float64[]

	# diagonal elements
	for (i, nd) in enumerate(node) 
		push!(data, nd.area*nd.Sy/dt + knobs.gamma*nd.sigma_w)		# water-water
		push!(data, knobs.gamma * alpha * nd.sigma_w) 				# water-NAPL
		push!(data, knobs.gamma * nd.sigma_n)	 					# NAPL-water	
		push!(data, nd.area*nd.Sy/dt + knobs.gamma*nd.sigma_n) 		# NAPL-NAPL
	end	
	
	# non-diagonal elements
	for (i, nd) in enumerate(node) 
		for cn in nd.connect_list
			push!(data, -knobs.gamma * connect[cn[1]].conduct_w) 			# water-water
			push!(data, -knobs.gamma * alpha * connect[cn[1]].conduct_w) 	# water-NAPL	
			push!(data, -knobs.gamma * connect[cn[1]].conduct_n) 			# NAPL-water	
			push!(data, -knobs.gamma * connect[cn[1]].conduct_n) 			# NAPL-NAPL
		end
	end
	
	return data
	
end


function RHSvector(connect::Array{Connect,1}, node::Array{Node,1}, num_nodes::Int64, alpha::Float64)

	# construct explicit matrix (run for each time step)
	b = zeros(Float64, 2*num_nodes)

	for i = 1:num_nodes
		b[i] = node[i].Qw - node[i].sigma_w * (node[i].z0 + node[i].hw + alpha*node[i].hn)
		b[i+num_nodes] = node[i].Qn - node[i].sigma_n * (node[i].z0 + node[i].hw + node[i].hn)		
	end	
	
	for cn in connect
	
		# water-balance
        b[cn.node_1] += cn.conduct_w *
			(node[cn.node_2].z0 + node[cn.node_2].hw + alpha*node[cn.node_2].hn)
        b[cn.node_2] += cn.conduct_w *
			(node[cn.node_1].z0 + node[cn.node_1].hw + alpha*node[cn.node_1].hn)

		# NAPL-balance
        b[cn.node_1 + num_nodes] += cn.conduct_n *
			(node[cn.node_2].z0 + node[cn.node_2].hw + node[cn.node_2].hn)
        b[cn.node_2 + num_nodes] += cn.conduct_n *
			(node[cn.node_1].z0 + node[cn.node_1].hw + node[cn.node_1].hn)
		
	end	
	
	return b

end

function MapNodeConnects(node::Array{Node, 1}, connect::Array{Connect, 1})
	# create list of node connections per each volume element; will be used to construct matrix of flow equations
	for (i, cn) in enumerate(connect)
		push!(node[cn.node_1].connect_list, (i, cn.node_2))
		push!(node[cn.node_2].connect_list, (i, cn.node_1))
	end
	return node
end

function UpdateTotConduct(node::Array{Node, 1}, connect::Array{Connect, 1})
	# update node total conductances (i.e., across all connections) for water and NAPL phases
	for nd in node
		nd.sigma_w = 0. 								# clear old values
		nd.sigma_n = 0.
	end
	for cn in connect
		node[cn.node_1].sigma_w += cn.conduct_w
		node[cn.node_2].sigma_w += cn.conduct_w
		node[cn.node_1].sigma_n += cn.conduct_n
		node[cn.node_2].sigma_n += cn.conduct_n		
	end
	return node
end


### main script ###


function napl_layer(t_end::Float64, mon_well::Int64)

	node, num_nodes = ReadNodes() 						# read nodes file and organize
	connect, num_connect = ReadConnects(node) 			# read connections file and organize
	node = MapNodeConnects(node, connect) 				# create lists of connected nodes, per volume element
	node = UpdateTotConduct(node, connect)				# update total conductances, per node
	fluid = ReadFluids() 								# read fluid properties
	alpha = fluid[2].rho/fluid[1].rho 					# fluid density ratio
	knobs = GetKnobs() 									# read in computational parameters
	
	t = 0.
	dt = knobs.dt_init
	monitor = Monitor[]
	
	# set up left-hand-side of flux balance equations matrix (conductances calc'd with initial conditions)
	data, row_index, col_index = LHSmatrix(connect, node, knobs, dt, num_nodes, alpha)
	
    while (t < t_end)

        complete = false

        while complete == false

            data = UpdateLHS(connect, node, knobs, dt, num_nodes, alpha) 		# update LHS elements (new dt, new conductances)
 			A = sparse(row_index, col_index, data, 2*num_nodes, 2*num_nodes)	# update sparse equation matrix
            b = RHSvector(connect, node, num_nodes, alpha)						# construct explicit vector
			global dh = \(A, b)	 												# solve equation set	
			
            # check maximum head change criterion (applied to both water and NAPL) at this time step size
			sum_complete = 0
			for i = 1:2*num_nodes
				sum_complete += (abs(dh[i]) > knobs.dh_max)
			end
			complete = 1 - sign(sum_complete)
            if complete == false 				# reduce time step size concentration change criterion not satisfied
				dt *= knobs.dt_decrease
				assert(dt > knobs.dt_min)
			end
			
		end	
			
        # update values
        t += dt 																# simulation time
		for (i, nd) in enumerate(node) 											# heads
			nd.hw += dh[i]
			nd.hn += dh[i+num_nodes]
		end
		push!(monitor, Monitor(t, node[mon_well].hw, node[mon_well].hn)) 		# monitoring well data logger
        
        # update time step
        dt *= knobs.dt_increase
        dt = min(dt, knobs.dt_max, t_end - t)	

		# update conductances (prior to next time step)
		for cn in connect
			cn.conduct_w, cn.conduct_n = Conduct(node, cn.node_1, cn.node_2, cn.dx, cn.len_inf)
		end
		
		# update total conductances, per node
		node = UpdateTotConduct(node, connect)
		
	end
	
	# write to model-wide output file
	fname = "final_state.csv"
	csvfile = open(fname,"w")
	line_out = "node" * "," * "x" * "," * "y" * "," * "hw" * "," * "hn"
	println(csvfile, line_out)	
	for i = 1:num_nodes
		line_out = string(i) * "," * string(node[i].x) * "," * string(node[i].y) * "," * string(node[i].hw) * "," * string(node[i].hn)
		println(csvfile,line_out)
	end
	close(csvfile)		
	println("Wrote final model state to output file.")

	# write to monitoring well file
	fname = "monitor.csv"
	csvfile = open(fname,"w")
	line_out = "t" * "," * "hw" * "," * "hn"
	println(csvfile, line_out)	
	for sample in monitor
		line_out = string(sample.time) * "," * string(sample.hw) * "," * string(sample.hn)
		println(csvfile,line_out)
	end
	close(csvfile)		
	println("Wrote monitoring well time series to output file.")	
	
	println("Done.")

end


### run script as napl_layer(t_end, mon_well)


# t_end = model end-time
# mon_well = node index number corresponding to monitor location

napl_layer(1.0, 1)
