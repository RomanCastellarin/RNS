#include <cstdio>
#define MOD(a,b) (((a)%(b)+(b))%(b))


class RNS{
    /// RESIDUAL NUMBER
    
    using ll = long long;

    public:
        // constructors
        RNS() : r{} {}
        RNS(ll x){
            for(int i = 0; i < M_SIZE; ++i)
                r[i] = MOD(x, M[i]);
        }
        // castings
        operator long long() const {
            ll S[M_SIZE] = { r[0] };
            ll a[M_SIZE] = { r[0] };
            for(int i = 1; i < M_SIZE; ++i){
                ll prod = 1;
                for(int j = 0; j < i; ++j)
                    prod = MOD(prod * powmod(M[j], M[i]-2, M[i]), M[i]);
                a[i] = MOD( MOD(r[i] - S[i-1], M[i]) * prod , M[i]);
                S[i] = S[i-1];
                if( a[i] ){
                    prod = 1;
                    for(int j = 0; j < i; ++j)
                        prod *= M[j];
                    S[i] += a[i] * prod;
                }
            }
            return S[M_SIZE-1];
        }
        
        // operators
        RNS operator + (RNS b) const {
            for(int i = 0; i < M_SIZE; ++i)
                b.r[i] = MOD(r[i] + b.r[i], M[i]);
            return b;
        }
        RNS operator - (RNS b) const {
            for(int i = 0; i < M_SIZE; ++i)
                b.r[i] = MOD(r[i] - b.r[i], M[i]);
            return b;
        }
        RNS operator * (RNS b) const {
            for(int i = 0; i < M_SIZE; ++i)
                b.r[i] = MOD(r[i] * b.r[i], M[i]);
            return b;
        }
        RNS operator / (RNS b) const { // fails when b and M[i] aren't coprime for some i. failure probability = 0.000000021
            for(int i = 0; i < M_SIZE; ++i)
                b.r[i] = MOD(r[i] * powmod(b.r[i], M[i]-2, M[i]), M[i]);
            return b;
        }

    private:
        // constants
        static const int M_SIZE = 6;
        static constexpr ll M[M_SIZE] = { 179426549, 181409623, 196659293, 254352227, 1779598819, 1779598963 };
        
        // helpers
        static ll powmod(ll b, ll e, ll mod) {
            if( e ){
                ll res = powmod(b, e>>1, mod);
                return (e & 1 ? (res*res) % mod * b : res*res) % mod;
            }else return 1;
        }

        // DATA
        ll r[M_SIZE];
};

constexpr RNS::ll RNS::M[];


int main(){
    RNS a = 768451336LL, b = 1028451213LL;
    RNS c = (a+b)*(a+b) + a*a*b*b*b - a*a - a*b - b*b ;
    RNS d = a*a*b*b*b;
    long long x = (d+d+d)/d; // = 3
    long long y = c-d; // = 790314708640670568 = ab

    
    printf("%Ld\n", x);
    printf("%Ld\n", y);

    return 0;
}
