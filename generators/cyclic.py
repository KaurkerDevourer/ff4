n = 7
mod = 1000000007
# generator for my lib


print("std::vector<Monomial<PrimeField<" + str(mod) + ">>> mon0;")
print()

kk = ""
for i in range(n):
    kk += "1"
    if i != n - 1:
        kk += ", "

print("mon0.push_back(Monomial(Term({" + kk + "}), PrimeField<" + str(mod) + ">(1)));")
print("mon0.push_back(Monomial(Term({0}), PrimeField<" + str(mod) + ">(-1)));")
print()
print("Polynomial<PrimeField<" + str(mod) + ">, GrevLexComp> p0(std::move(mon0));")
print()

for i in range(1, n):
    print("std::vector<Monomial<PrimeField<" + str(mod) + ">>> mon" + str(i) + ";")
    print()
    l = i
    kkk = []
    arrr = []
    for k in range(n): # сдвиг
        kk = ""
        arr = [0] * n
        for t in range(l): # проставляем все единички
            arr[(k + t) % n] = 1
        arrr.append(arr[::-1])

    arrr.sort()

    for arr in arrr:
        arr = arr[::-1]
        while arr[len(arr)- 1] == 0:
            arr.pop()

        kk = ""
        for t in range(len(arr)): # записываем в строку
            kk += str(arr[t])
            if t != len(arr) - 1:
                kk += ", "
        
        print("mon" + str(i) + ".push_back(Monomial(Term({" + kk + "}), PrimeField<" + str(mod) + ">(1)));")

    print()
    print("Polynomial<PrimeField<" + str(mod) + ">, GrevLexComp> p" + str(i) + "(std::move(mon" + str(i) + "));")
    print()

finish = ""
for i in range(n):
    finish += "p" + str(i)
    if i != n - 1:
        finish += ", "
print("TPolynomials<PrimeField<" + str(mod) + ">, GrevLexComp> test = {" + finish + "};")
