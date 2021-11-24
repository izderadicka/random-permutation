use std::fmt::{Debug, Display};

use anyhow::Result;
use glass_pumpkin::prime;
use log::debug;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{Num, NumAssignOps, NumOps, One};
use rand::{Rng, SeedableRng};

pub struct Generator<T> {
    limit: T,
    modulo: T,
    multiplier: T,
    state: T,
    count: T,
}

enum RandGen {
    Seeded(rand::rngs::StdRng),
    Unseeded(rand::rngs::ThreadRng),
}

impl RandGen {
    fn new(seed: Option<impl AsRef<[u8]>>) -> Self {
        match seed {
            Some(s) => {
                let mut seed = [0u8; 32];
                let bytes = s.as_ref();
                let sz = seed.len().min(bytes.len());
                seed[..sz].copy_from_slice(&bytes[..sz]);
                RandGen::Seeded(rand::rngs::StdRng::from_seed(seed))
            }
            None => RandGen::Unseeded(rand::thread_rng()),
        }
    }

    fn gen_range<T, R>(&mut self, range: R) -> T
    where
        T: rand::distributions::uniform::SampleUniform,
        R: rand::distributions::uniform::SampleRange<T>,
    {
        match self {
            RandGen::Seeded(r) => r.gen_range(range),
            RandGen::Unseeded(r) => r.gen_range(range),
        }
    }
}

/// Simple generation of permutation in memory using Fisher-Yates shuffles (also known as Knuth algorithm)
/// generates numbers 1 to limit (inclusive)
pub fn random_permutation<T, S>(limit: T, seed: Option<S>) -> Vec<T>
where
    T: Integer + Clone,
    usize: From<T>,
    S: AsRef<[u8]>,
{
    let size: usize = limit.into();
    let mut p: Vec<T> = Vec::with_capacity(size);
    (1..=size).fold(T::one(), |item, _| {
        p.push(item.clone());
        item + T::one()
    });
    let mut rng = RandGen::new(seed);
    for i in 0..(size - 1) {
        let j = i + rng.gen_range(0..size - i);
        p.swap(i, j)
    }
    p
}

impl<T> Generator<T>
where
    for<'a> T: Integer + Clone + NumOps<&'a T> + NumAssignOps<T>,
{
    pub fn new(modulo: T, multiplier: T) -> Self {
        Self::with_limit_and_state(modulo.clone(), multiplier, T::one(), modulo)
    }

    pub fn with_limit_and_state(modulo: T, multiplier: T, state: T, limit: T) -> Self {
        assert!(state > T::zero() && state < modulo);
        assert!(limit > T::one() && limit <= modulo);
        Self {
            limit,
            modulo,
            multiplier,
            state,
            count: T::one(),
        }
    }

    fn step(&mut self) {
        self.state = self.multiplier.clone() * &self.state % &self.modulo;
        self.count += T::one();
    }
}

impl<T> Iterator for Generator<T>
where
    for<'a> T: Integer + Clone + NumOps<&'a T> + NumAssignOps<T>,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        while self.count < self.modulo {
            self.step();
            if self.state < self.limit {
                return Some(self.state.clone());
            }
        }

        None
    }
}

/// Finds first prime after given limit
pub fn generate_prime<T>(limit: T) -> Result<T>
where
    T: Into<BigUint> + Num,
{
    let limit: BigUint = limit.into();
    let res = generate_prime_int(limit);
    T::from_str_radix(&res.to_str_radix(10), 10)
        .map_err(|_| anyhow::Error::msg("Invalid result, probably too big for given type"))
}

/// Finds first prime after given limit
fn generate_prime_int(limit: impl Into<BigUint>) -> BigUint {
    let mut n = limit.into();
    let two = BigUint::from(2u8);
    if n.is_multiple_of(&two) {
        n += BigUint::one();
    } else {
        n = n + &two;
    }

    loop {
        if prime::check(&n) && prime::strong_check(&n) {
            return n;
        }
        n += &two;
    }
}

fn powmod<N>(mut x: N, mut y: N, p: &N) -> N
where
    for<'a> N: Integer + Clone + NumOps<&'a N> + NumAssignOps<&'a N> + NumAssignOps<N>,
{
    if y.is_zero() {
        N::one()
    } else {
        let mut res = N::one();
        x %= p;
        let two = N::one() + N::one();
        while y > N::zero() {
            if y.is_odd() {
                res = res.clone() * &x % p;
                y -= N::one();
            } else {
                x = x.clone() * &x % p;
                y /= &two;
            }
        }

        res
    }
}

/// finds bigest primitive root for given n
/// n must be prime bigger then 7 - small primes are not interesting for this use case
pub fn find_primitive_root<N>(p: N) -> Option<N>
where
    for<'a> N:
        Integer + Clone + NumOps<&'a N> + NumAssignOps<&'a N> + NumAssignOps<N> + Display + Debug,
{
    let phi = p.clone() - N::one();
    let mut n = phi.clone();

    // Find prime factors
    let mut divs = Vec::new();
    let mut i = N::one() + N::one();
    let giga = N::from_str_radix("1000000000", 10).ok();
    while i.clone() * &i <= n {
        if n.is_multiple_of(&i) {
            debug!("Found factor : {}", i);
            divs.push(i.clone());
            n /= &i;
            while n.is_multiple_of(&i) {
                n /= &i;
            }
        }
        if let Some(ref giga) = giga {
            if i.is_multiple_of(giga) {
                debug!("Processed 1G, remains {}", n.clone() / giga.clone());
            }
        }
        i += N::one();
    }
    if n > N::one() {
        divs.push(n)
    }

    debug!("Found  {} prime factors : {:?}", divs.len(), divs);

    let two = N::one() + N::one();
    let max = phi.clone() - N::one();
    let mut res = max / &two; // start app from middle
    while res >= two {
        let mut ok: bool = true;
        for i in &divs {
            let rem = powmod(res.clone(), phi.clone() / i, &p);
            if rem == N::one() {
                ok = false;
                break;
            }
        }
        if ok {
            debug!("Prim. root: {}", res);
            //roots.push(res.clone());
            return Some(res);
        }
        res -= N::one();
    }
    //debug!("All roots: {:?}", roots);
    //roots.get(0).cloned()
    None
}

pub fn pi<T: Integer>(s: impl AsRef<str>) -> Result<T> {
    T::from_str_radix(s.as_ref(), 10).map_err(|_| anyhow::Error::msg("Invalid lower_limit"))
}

#[cfg(test)]
mod tests {
    use std::{collections::HashSet, fmt::Display, hash::Hash};

    use num_traits::Num;

    use super::*;

    fn is_root<N>(root: N, p: N, expected_size: usize) -> bool
    where
        for<'a> N: Integer
            + Clone
            + NumOps<&'a N>
            + NumAssignOps<&'a N>
            + NumAssignOps<N>
            + Hash
            + Display,
    {
        let mut res = HashSet::new();
        let mut x = N::one();
        while x < p {
            let y = powmod(root.clone(), x.clone(), &p);
            res.insert(y);
            x += N::one();
        }
        if res.len() == expected_size {
            true
        } else {
            false
        }
    }

    #[test]
    fn test_powmod() {
        let res: Vec<_> = (1..=6u32).map(|x| 5u32.pow(x) % 7).collect();
        for i in 1..=6u32 {
            let x = powmod(5u32, i, &7);
            assert_eq!(res[(i - 1) as usize], x)
            //println!("{} : {}<=>{}", i, res[(i-1) as usize], x);
        }
    }

    #[test]
    fn test_prim_root() {
        let n1 = 13u8;
        let r1 = find_primitive_root(n1).expect("root exists");
        println!("Primitive root for {} is {}", n1, r1);
        assert!(is_root(r1, n1, (n1 - 1) as usize));
        let n2_num = 104729usize;
        let n2 = BigUint::from(n2_num);
        let r2 = find_primitive_root(n2.clone()).expect("root exists");
        println!("Primitive root for {} is {}", n2, r2);
        assert!(is_root(r2.clone(), n2.clone(), n2_num - 1));
    }
    #[test]
    fn test_prime_m607() {
        let prime="531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127";
        let prime = BigUint::from_str_radix(prime, 10).unwrap();
        assert!(prime::check(&prime));
        assert!(prime::strong_check(&prime))
    }

    #[test]
    fn test_generated_prime() {
        let prime = generate_prime_int(1_000_000u32);
        println!("{}", prime);
        assert!(prime::check(&prime));
        assert!(prime::strong_check(&prime))
    }

    #[test]
    fn test_generator() {
        let fg = Generator::new(13u32, 11);
        let nums: Vec<_> = fg.collect();
        assert_eq!(12, nums.len());
        let nums: HashSet<_> = nums.into_iter().collect();
        assert_eq!(12, nums.len());
        assert!(nums.iter().all(|x| *x >= 1 && *x <= 12))
    }

    #[test]
    fn test_limited_generator() {
        let fg = Generator::with_limit_and_state(10007u32, 10005, 500, 5001);
        let nums: HashSet<_> = fg.collect();
        assert_eq!(5000, nums.len());
        assert!(nums.iter().all(|x| *x >= 1 && *x <= 5000))
    }

    #[test]
    fn test_bigint_generator() {
        let p = pi("7138297498312732814783179").unwrap();
        assert!(prime::strong_check(&p));
        // TODO Why this number is not working if this is primitive root by definition
        //let g = bi("7138297498312732814783178").unwrap();
        let g: BigUint = pi("3569148749156366407391588").unwrap();
        println!("diff {}", p.clone() - g.clone());
        let fg = Generator::with_limit_and_state(p.clone(), g, pi("999").unwrap(), p);
        let nums: HashSet<_> = fg.take(100000).collect();
        assert_eq!(100000, nums.len());
    }

    #[ignore]
    #[test]
    fn test_bigint_generator2() {
        let p = pi("71382991").unwrap();
        assert!(prime::strong_check(&p));
        let g: BigUint = pi("71382986").unwrap();
        println!("diff {}", p.clone() - g.clone());
        let fg =
            Generator::with_limit_and_state(p.clone(), g, pi("999").unwrap(), pi("1000").unwrap());
        let nums: HashSet<_> = fg.collect();
        assert_eq!(999, nums.len());
    }

    #[test]
    fn test_random_permutation() {
        let p = random_permutation::<_, &str>(100u8, None);
        let nums: HashSet<_> = p.into_iter().collect();
        assert_eq!(100, nums.len());
        assert!(nums.iter().all(|x| *x >= 1 && *x <= 100))
    }
}
