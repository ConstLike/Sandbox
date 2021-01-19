"""  Ratio(n, d)
return fraction with nomerator (n) and denominator (d) """
struct Ratio <: Real
	numer::Int
	denom::Int

	function Ratio(n, d)
		let s = (d >= 0 ? 1 : -1), g = gcd(n, d)
			return new(s * n ÷ g, abs(d) ÷ g)
		end
	end
end

begin
    numer(r::Ratio) = r.numer
    denom(r::Ratio) = r.denom

    function Base.:+(x::Ratio, y::Ratio)
    	Ratio(numer(x) * denom(y) + numer(y) * denom(x),
              denom(x) * denom(y))
    end
    function Base.:+(x::Ratio, y::Integer)
        let y = Ratio(y, 1)
            Ratio(numer(x) * denom(y) + numer(y) * denom(x),
                  denom(x) * denom(y))
        end
    end
    function Base.:+(x::Integer, y::Ratio)
        let x = Ratio(x, 1)
            Ratio(numer(x) * denom(y) + numer(y) * denom(x),
                  denom(x) * denom(y))
        end
    end

    function Base.:-(x::Ratio, y::Ratio)
        Ratio(numer(x) * denom(y) - numer(y) * denom(x),
              denom(x) * denom(y))
    end
    function Base.:-(x::Ratio, y::Integer)
        let y = Ratio(y, 1)
            Ratio(numer(x) * denom(y) - numer(y) * denom(x),
                  denom(x) * denom(y))
        end
    end
    function Base.:-(x::Integer, y::Ratio)
        let x = Ratio(x, 1)
            Ratio(numer(x) * denom(y) - numer(y) * denom(x),
                  denom(x) * denom(y))
        end
    end

    function Base.:*(x::Ratio, y::Ratio)
    	Ratio(numer(x) * numer(y),
              denom(x) * denom(y))
    end
    function Base.:*(x::Integer, y::Ratio)
        Ratio(x * numer(y),
              denom(y))
    end
    function Base.:*(x::Ratio, y::Integer)
        Ratio(numer(x) * y,
              denom(x))
    end

    function Base.:/(x::Ratio, y::Ratio)
    	Ratio(numer(x) * denom(y),
              denom(x) * numer(y))
    end
    function Base.:/(x::Integer, y::Ratio)
    	Ratio(x * denom(y),
              numer(y))
    end
    function Base.:/(x::Ratio, y::Integer)
        Ratio(numer(x),
              denom(x) * y)
    end

    function Base.:(==)(x::Ratio, y::Ratio)
    	numer(x) * denom(y) == numer(y) * denom(x)
    end
    function Base.:(==)(x::Ratio, y::Integer)
        let y = Ratio(y, 1)
            numer(x) * denom(y) == numer(y) * denom(x)
        end
    end
    function Base.:(==)(x::Integer, y::Ratio)
        let x = Ratio(x, 1)
            numer(x) * denom(y) == numer(y) * denom(x)
        end
    end

    function Base.:abs(x::Ratio)
        Ratio(abs(numer(x)), abs(denom(x)))
    end

    Base.show(io::IO, r::Ratio) = print(io, numer(r), '/', denom(r))
end

Ratio(1, 2) + Ratio(3, 4)
Ratio(1, -2)
Ratio(5, -2) + Ratio(1, -2)
