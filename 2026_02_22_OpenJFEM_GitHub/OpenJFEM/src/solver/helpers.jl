# helpers.jl — Small utility functions for the solver

function log_msg(message::String)
    println("[$(Dates.format(now(), "HH:MM:SS"))] $message")
end

# CBAR V-vector resolution (G0 grid point or direct vector)
function resolve_bar_vref(bar, p_ga::SVector{3,Float64}, id_map, node_coords)
    g0 = get(bar, "G0", 0)
    if g0 > 0 && haskey(id_map, g0)
        ig0 = id_map[g0]
        p_g0 = SVector{3}(node_coords[ig0,1], node_coords[ig0,2], node_coords[ig0,3])
        return p_g0 - p_ga
    else
        v_raw = bar["V"]
        return SVector{3}(v_raw...)
    end
end

# Skew-symmetric cross-product matrix: skew(w) * v = w × v
@inline function skew3(w)
    return @SMatrix [0.0 -w[3] w[2]; w[3] 0.0 -w[1]; -w[2] w[1] 0.0]
end

# Get CBAR offset vectors and effective endpoints
function bar_offsets_and_endpoints(bar, p1::SVector{3,Float64}, p2::SVector{3,Float64})
    wa_raw = get(bar, "WA", [0.0, 0.0, 0.0])
    wb_raw = get(bar, "WB", [0.0, 0.0, 0.0])
    wa = SVector{3}(wa_raw...)
    wb = SVector{3}(wb_raw...)
    has_offset = (wa[1]^2+wa[2]^2+wa[3]^2) > 1e-20 || (wb[1]^2+wb[2]^2+wb[3]^2) > 1e-20
    p1_eff = has_offset ? p1 + wa : p1
    p2_eff = has_offset ? p2 + wb : p2
    return wa, wb, has_offset, p1_eff, p2_eff
end

# Fast shell element frame computation
function shell_element_frame_fast(p1::SVector{3,Float64}, p2::SVector{3,Float64}, p3::SVector{3,Float64}, p4::SVector{3,Float64}, n::Int)
    if n == 4
        d13 = p3 - p1; d24 = p4 - p2
        v3 = normalize(cross(d13, d24))
        d13n = normalize(d13); d24n = normalize(d24)
        x_raw = d13n - d24n
        if dot(x_raw, x_raw) < 1e-20; x_raw = d13n + d24n; end
        x_raw = normalize(x_raw)
        v1 = normalize(x_raw - dot(x_raw, v3)*v3)
        v2 = cross(v3, v1)
        return v1, v2, v3
    else
        v1 = normalize(p2 - p1)
        v3 = normalize(cross(v1, p3 - p1))
        v2 = cross(v3, v1)
        return v1, v2, v3
    end
end

function rotation_from_normal(n_avg::Vector{Float64})
    z = normalize(n_avg)
    ref = abs(z[3]) < 0.9 ? SVector(0.0, 0.0, 1.0) : SVector(1.0, 0.0, 0.0)
    x = normalize(cross(ref, SVector(z...)))
    y = cross(SVector(z...), x)
    return hcat(x, y, z)
end

function get_coord_transform(model, cid, vec)
    if cid == 0; return vec; end
    if !haskey(model["CORDs"], string(cid)); return vec; end
    cord = model["CORDs"][string(cid)]
    R = hcat(cord["U"], cord["V"], cord["W"])
    return R * vec
end
