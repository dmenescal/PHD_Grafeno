##################################################################################################
# Execução principal da extração completa de dados (Cones de luz com e sem nanofitas de grafeno) #
##################################################################################################

############################################
# Parâmetros iniciais: Indices de Refração:
############################################
n1 = 1;
n2 = (3.42)^2;
n3 = (1.98)^2;

# Variação de cada valor imaginário para n2
delta = 0.1;
# Vetor de valores para a parte imaginária para n2
n2_irange = [0.1:delta:10*n2];

# For testing, uncomment
extractData(n1, n2+9.8*i, n3);

####################################
# Extração das Bases sem nanofitas:
####################################

for n = 1:length(n2_irange)
  #n2_complex = n2 + n2_irange(n)*i;
  #extractData(n1, n2_complex, n3);
end

####################################
# Extração das Bases com nanofitas:
####################################

for n = 1:length(n2_irange)
  #n2_complex = n2 + n2_irange(n)*i;
  #extractGraphData(n2_complex, n3);
end
