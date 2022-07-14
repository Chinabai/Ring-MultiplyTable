//
//  main.cpp
//  Ring_Multiply
//
//  Created by Chinabai on 2022/7/14.
//

#include <iostream>
#include <algorithm>
#include <cmath>
#include <bitset>
#include <fstream>
#include <string>
#include <sys/uio.h>
#include <stdio.h>
#include <vector>
#include <sys/time.h>
using namespace std;

using LL = long long;

/**
 多项式类
 */
class Polynomial {
    public:
        LL data;
    /**
     获取二进制的位数
     */
    const size_t binary_length()const {
        LL tmp = data;
        int len = 0;
        while(tmp) {
            len++;
            tmp >>= 1;
        }
        return len;
    }
    
    /**
     构造函数
     */
    Polynomial(LL data) {
        this->data = data;
    }
    
    /**
     模2加法
     */
    Polynomial operator + (const Polynomial &o) {
        return data ^ o.data;
    }
    
    /**
     模2减法
     */
    Polynomial operator - (const Polynomial &o) {
        return data ^ o.data;
    }
    
    /**
     模2乘法
     */
    Polynomial operator * (const Polynomial &o) {
        Polynomial ans = Polynomial(0);
        Polynomial tmp = o;
        for(int i = 0; tmp.data; i++, tmp.data >>= 1) {
            if(tmp.data & 1) {
                ans = ans + Polynomial(this->data << i);
            }
        }
        return ans;
    }
    
    /**
     模2除法
     */
    Polynomial operator / (const Polynomial &o) {
        size_t len = this->binary_length();
        size_t o_len = o.binary_length();
        int width = (int)len - (int)o_len;
        Polynomial remainder = *this;
        Polynomial quotient = Polynomial(0);
        while(width >= 0) {
            if(remainder.data >> (o_len + width - 1)) {
                quotient.data = quotient.data + (1LL << width);
                remainder = remainder - Polynomial(o.data << width);
            }
            width--;
        }
        return quotient;
    }
    
    /**
     模2模
     */
    Polynomial operator % (const Polynomial &o) {
        size_t len = this->binary_length();
        size_t o_len = o.binary_length();
        int width = (int)len - (int)o_len;
        Polynomial remainder = *this;
        Polynomial quotient = Polynomial(0);
        while (width >= 0) {
            if(remainder.data >> (o_len + width - 1)) {
                quotient.data = quotient.data + (1LL << width);
                remainder = remainder - Polynomial(o.data << width);
            }
            width--;
        }
        return remainder;
    }
    
    /**
     重载输出
     */
    friend ostream &operator << (ostream &os, const Polynomial &o);
};

/**
 重载输出运算符
 */
ostream &operator << (ostream &os, const Polynomial &o) {
    auto data = o.data;
    if(data == 0) {
        os << 0;
        return os;
    }
    bool flag = false;
    for(int i = 0; data; i++, data >>= 1) {
        if(data & 1) {
            if(flag) {
                os << "+";
            }
            flag = true;
            //i = 0 时 为 1
            if(i) {
                os << "x^" << i;
            } else {
                os << 1;
            }
        }
    }
    return os;
}

/**
 欧几里德求逆元
 */
pair<Polynomial, Polynomial> exgcd(Polynomial a, Polynomial b) {
    if(b.data == 0) {
        return {Polynomial(1), Polynomial(0)};
    }
    auto factor = exgcd(b, a%b);
    auto &x = factor.first;
    auto &y = factor.second;
    auto xp = y, yp = x - (a / b * y);
//    cout << "x = " << x << ", ";
//    cout << "y = " << y << endl;
    return {xp, yp};
}

/**
 tools
 */

/**
 十进制和十六进制转换
 */
string getHex(Polynomial b) {
    LL x, deci = b.data;
    string result = "";
    char mid;
    if(deci == 0) {
        return "0";
    }
    while(deci != 0) {
        x = deci % 16;
        if(x < 10) {
            mid = x + '0';
        } else {
            mid = x + 'A' - 10;
        }
        result = mid + result;
        deci = deci / 16;
    }
    return result;
}

/**
 求素数
 求解n以内(含n)的素数
 */
vector<LL> getPrime(LL n) {
    //标记数组,flag[i]==0表示i为素数,flag[i]==1表示i为合数
    bool flag[n + 1];
    memset(flag, 0, sizeof(flag));
    // 素数个数
    vector<LL> prime;
    int cnt = 0;
    for (int i = 2; i <= n; ++i) {
        if (!flag[i]) {
            //将i加入素数表
            prime.push_back(i);
            cnt++;
        }
        for (int j = 0; j < cnt; ++j){
            //保证每个合数只会被它的最小质因数筛去
            if (i * prime[j] > n)  break;
            flag[i * prime[j]] = 1;
            if (i % prime[j] == 0)  break;
        }
    }
    return prime;
}

/**
 求乘法表
 在模m(x)上
 */
Polynomial getMultiply(LL a, LL b, Polynomial M) {
    Polynomial mul = Polynomial(a) * Polynomial(b);
    return mul % M;
}

/**
 生成p = M的乘法表
 在模m(x)上
 */
void getMultiplyTable(Polynomial M) {
    LL bit_size = M.binary_length();
    LL table_size = (int)M.data / 16 + 1;
    //创建表
    ofstream ofile;
    ofile.open("multiply_p"+to_string(bit_size)+".csv", ios::out | ios::trunc);
    //设置表头
    ofile << M << ",,";
    for(LL i = 1; i < table_size * 16; i++) {
        string table_head = getHex(Polynomial(i));
        if(i % 16 == 0 && i != 0) {
            ofile << ",,";
           // cout << endl;
        }
        ofile << table_head << ",";
    }
  
    //表内容
    for(LL i = 1; i < M.data; i++) {
        ofile << "\n," << getHex(Polynomial(i)) << ",";
        for(LL j = 1; j < M.data; j++) {
            string mul_result = getHex(getMultiply(i, j, M));
            if(j % 16 == 0) {
                ofile << ",,";
            }
            ofile << mul_result << ",";
        }
    }
    ofile.close();
}

/**
 生成p = 5~M
 的乘法表
 在模m(x)上
 */
void getMultiplyTables(Polynomial M) {
    vector<LL> primes = getPrime(M.data);
    Polynomial o = Polynomial(0);
    for(int i = 2; i < primes.size(); i++) {
        o.data = (1 << primes[i]) - 1 ;
        getMultiplyTable(o);
    }
}


int main(int argc, const char * argv[]) {
    LL p = 7;
    //生成乘法表 p = 5 ～ p范围内的
    getMultiplyTables(Polynomial(p));
    return 0;
}
