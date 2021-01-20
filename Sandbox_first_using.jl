"""  Ratio(n, d)
return fraction with nomerator (n) and denominator (d) """
struct Ratio <: Real
	numer::Int
	denom::Int

	function Ratio(n, d)
		let s = (d >= 0 ? 1 : -1), g = gcd(n, d)
            if d == 0
                return error("Argument of 'denominator' cannot be zero")
            else
                return new(s * n ÷ g, abs(d) ÷ g)
            end
		end
	end
end

"""  Mixed_Frac(r, k)
ruturn an integer part 'r' of the improper fraction 'r' """
struct Mixed_Frac <: Real
    num::Int
	den::Int
    ing::Int

	function Mixed_Frac(r::Ratio, itg=0)
        let n = r.numer, d = r.denom,
			s = (d >= 0 ? (n >= 0 ? 1 : -1) : (n >= 0 ? -1 : 1)),
            d = abs(d), itg = abs(itg), n = abs(n), g = div(n, d), itg = s * g
            if d > n
                return new(s * n, d, itg)
            else
                return new(s * (n - g * d), d, itg)
            end
        end
    end
end


begin
    numer(r::Ratio) = r.numer
    denom(r::Ratio) = r.denom
    num(r::Mixed_Frac) = r.num
    den(r::Mixed_Frac) = r.den
    ing(r::Mixed_Frac) = r.ing

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
    function Base.:^(x::Ratio, y::Integer)
        Ratio(numer(x)^y,
              denom(x)^y)
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
    Base.show(io::IO, r::Mixed_Frac) =
        (ing(r) == 0 ?
            (num(r) == 0 ?
                print(io, num(r))
                : print(io, num(r), '/', den(r))
            )
            : (num(r) == 0 ?
                print(io, ing(r))
                : print(io, ing(r), " and ", num(r), '/', den(r))
            )
        )
end

""" Testing:         """
one_half = Ratio(1, 0) # WRONG: Argument of 'denominator' cannot be zero
one_half = Ratio(1, 2) # = 1/2
one_third = Ratio(1, 3) # = 1/3
minus_quarter = Ratio(1, -4) # = -1/4

one_half / minus_quarter # = -2/1
Mixed_Frac(one_half / minus_quarter) # = -2

Mixed_Frac(one_half * one_third / minus_quarter) # = -2/3
Mixed_Frac(one_half - one_third / minus_quarter) # = 1 and 5/6
Mixed_Frac(one_half * one_third * minus_quarter) # = -1/24
Mixed_Frac(48 * one_half * one_third * minus_quarter) # = -2
Mixed_Frac(3 * one_half * one_third) # = 1/2
Mixed_Frac(3 / one_half * one_third) # = 2
Mixed_Frac(3 + one_half * one_third) # 3 and 1/6
Mixed_Frac(one_half - one_half) # 0
Mixed_Frac(one_third ^ 3) # = 1/27
Mixed_Frac(abs(minus_quarter)) # = 1/4
