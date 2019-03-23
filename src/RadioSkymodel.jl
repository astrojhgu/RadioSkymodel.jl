module RadioSkymodel

import AstroLib
import Healpix

struct SkyModelData
    i_template::Array{Float64, 1}
    u_template::Array{Float64, 1}
    q_template::Array{Float64, 1}
    beta_template::Array{Float64, 1}
    nu_0_i::Float64
    nu_0_q::Float64
    nu_0_u::Float64
    resolution::Healpix.Resolution
end


function SkyModelData(i::Healpix.Map, q::Healpix.Map, u::Healpix.Map, beta::Healpix.Map, nu_0_i::Float64, nu_0_p::Float64)::SkyModelData
    SkyModelData(i.pixels, u.pixels, q.pixels, beta.pixels, nu_0_i, nu_0_p, nu_0_p, i.resolution)
end


function radec2gal(ra, dec)
    result = Tuple{Float64, Float64}[]
    for (r,d) in zip(ra, dec)
        push!(result, AstroLib.euler(r, d, 1; FK4 = false, radians = true))
    end
    collect.(zip(result...))
end

function radec2θφ(ra, dec)
    φ, b=radec2gal(ra, dec)
    θ = π / 2 .-b
    θ,φ
end

function get_iqu(ra::Array, dec::Array, freq::Array, sm::SkyModelData)::Tuple{Array, Array, Array}
    θ, φ=radec2θφ(ra, dec)
    out_i = zeros(Float64, size(ra, 0 + 1),size(freq, 0 + 1))
    out_q = zeros(Float64, size(ra, 0 + 1),size(freq, 0 + 1))
    out_u = zeros(Float64, size(ra, 0 + 1),size(freq, 0 + 1))
    for j in eachindex(ra)
        println(j, " ", size(ra,1))
        idx, wgt=Healpix.getinterpolRing(sm.resolution, θ[j], φ[j])
        i0=sum(sm.i_template[i]*w for (i,w) in zip(idx, wgt))
        q0=sum(sm.q_template[i]*w for (i,w) in zip(idx, wgt))
        u0=sum(sm.u_template[i]*w for (i,w) in zip(idx, wgt))
        beta=sum(sm.beta_template[i]*w for (i,w) in zip(idx, wgt))
        i=i0.*(freq./sm.nu_0_i).^beta
        q=q0.*(freq./sm.nu_0_q).^beta
        u=u0.*(freq./sm.nu_0_u).^beta
        out_i[j, :]=i
        out_q[j, :]=q
        out_u[j, :]=u
    end
    (out_i/1e6, out_q/1e6, out_u/1e6)
end

end # module
