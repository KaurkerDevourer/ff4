n = 9
# generator for my lib


print("std::vector<Monomial<PrimeField<31>>> mon0;")

kk = ""
for i in range(n):
    kk += "1"
    if i != n - 1:
        kk += ", "

print("mon0.push_back(Monomial(Term({" + kk + "}), PrimeField<31>(1)));")
print("mon0.push_back(Monomial(Term({0}), PrimeField<31>(-1)));")
print("Polynomial<PrimeField<31>, GrevLexComp> p0(std::move(mon0));")
print()

for i in range(1, n - 1):
    print("std::vector<Monomial<PrimeField<31>>> mon" + str(i) + ";")
    print()
    l = i + 1
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
        
        print("mon" + str(i) + ".push_back(Monomial(Term({" + kk + "}), PrimeField<31>(1)));")

    print()
    print("Polynomial<PrimeField<31>, GrevLexComp> p" + str(i) + "(std::move(mon" + str(i) + "));")
    print()

finish = ""
for i in range(n - 1):
    finish += "p" + str(i)
    if i != n - 2:
        finish += ", "
print("TPolynomials<PrimeField<31>, GrevLexComp> test = {" + finish + "};")
