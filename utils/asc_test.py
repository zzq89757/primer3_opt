
def asc_test() -> list:
    base_li = ['A','C','G','T','I','O','-']
    ord_li = [ord(x) for x in base_li]
    print(ord_li)
    mu_li = []
    for i in range(len(base_li)):
        for j in range(i,len(base_li)):
            print(base_li[i],end="->")
            print(base_li[j],end="\t")
            print(ord_li[i] * ord_li[j])
            if ord_li[i] * ord_li[j] in mu_li:
                print(base_li[i],end="*")
                print(base_li[j],end="\t")
                print(ord_li[i] * ord_li[j])
            mu_li.append(ord_li[i] * ord_li[j])
            
            # 整除I（73）
            if (ord_li[i] * ord_li[j] % 73 == 0):
                print(base_li[i],end="*")
                print(base_li[j],end="\t")
                print(ord_li[i] * ord_li[j])
            # print(ord_li[i] * ord_li[j] % 65)
            # print(ord_li[i] * ord_li[j] % 45)
            
            # AT 5460, CG 4757, others % 73 == 0 and % 45 != 0 (remove I-)


if __name__ == "__main__":
    asc_test()