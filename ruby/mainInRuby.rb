require 'bigdecimal'
# require 'Math'
BigDecimal.limit(10)
L = BigDecimal.new("1")
N = BigDecimal.new("10")
dx = BigDecimal.new("0")
dx = L/N # => #<BigDecimal:7faa3b890c78,'0.1E0',9(36)>
x = BigDecimal.new("0")
h = []
H = BigDecimal.new("1")
w = BigDecimal.new("1")
a = BigDecimal.new("0.2") # => #<BigDecimal:7faa3b890570,'0.2E0',9(18)>
a = a*H # => #<BigDecimal:7faa3b8900c0,'0.2E0',9(36)>
w = 0.1*L # => 0.1


## Loop for initial height
1.upto(N) do |grid_numbers|
  grid_numbers # => 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  x = (grid_numbers-0.5)*dx # => 0.05, 0.15000000000000002, 0.25, 0.35000000000000003, 0.45, 0.55, 0.65, 0.75, 0.8500000000000001, 0.9500000000000001
  h.push( H + a*Math.exp(-(x-L/2)**2/(w)**2) ) # h(x,0)
end







h # => [#<BigDecimal:7faa3b89a7f0,'0.1E1',9(36)>, #<BigDecimal:7faa3b8998a0,'0.1000000957E1',18(36)>, #<BigDecimal:7faa3b898978,'0.1000386091E1',18(36)>, #<BigDecimal:7faa3b8a3c10,'0.1021079845E1',18(36)>, #<BigDecimal:7faa3b8a2ec8,'0.1155760157E1',18(36)>, #<BigDecimal:7faa3b8a2180,'0.1155760157E1',18(36)>, #<BigDecimal:7faa3b8a1438,'0.1021079845E1',18(36)>, #<BigDecimal:7faa3b8a0510,'0.1000386091E1',18(36)>, #<BigDecimal:7faa3b8ab5f0,'0.1000000957E1',18(36)>, #<BigDecimal:7faa3b8aa6a0,'0.1E1',9(36)>]
