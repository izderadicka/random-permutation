use std::process::exit;

use anyhow::{bail, Error, Result};
use log::debug;
use num_bigint::BigUint;
use num_traits::One;
use rand::Rng;
use random_permutation::*;
use sha2::Digest;
use structopt::StructOpt;

#[derive(StructOpt)]
enum Actions {
    #[structopt(
        help = "Generate params for Lehman LCG - prime for modulo and primitive root for multiplier"
    )]
    GenerateParams {
        #[structopt(required = true, help = "Prime will be bigger then")]
        lower_limit: String,
        #[structopt(
            short = "b",
            long,
            help = "Use bigint - no upper limit, but slower (app. 7x)"
        )]
        use_bigint: bool,
    },
    #[structopt(help = "Prints psedorandom permutation")]
    Permutation {
        #[structopt(
            required = true,
            help = "Multiplier - must be primitive root of modulo parameter"
        )]
        multiplier: String,
        #[structopt(required = true, help = "Modulo - must be prime")]
        modulo: String,
        #[structopt(short, long, help = "Seed value, must be smaller then modulo")]
        seed: Option<String>,
        #[structopt(
            short,
            long,
            help = "Upper limit for permutation, must be smaller then modulo"
        )]
        upper_limit: Option<String>,
        #[structopt(short, long, help = "Print only first x numbers")]
        take: Option<usize>,
        #[structopt(
            short = "b",
            long,
            help = "Use bigint - no upper limit, but slower (app. 2x)"
        )]
        use_bigint: bool,
    },

    #[structopt(help = "Prints random using Fisher-Yates and CSPRG")]
    RandomPermutation {
        #[structopt(
            required = true,
            help = "Max element in permutation - permutating 1..=this"
        )]
        upper_limit: usize,
        #[structopt(short, long, help = "Seed value, any string value is ok - It'll be hasheded")]
        seed: Option<String>,
        #[structopt(short, long, help = "Print only first x numbers")]
        take: Option<usize>,
    },
}

fn flip<T, E>(x: Option<Result<T, E>>) -> Result<Option<T>, E> {
    x.map_or(Ok(None), |v| v.map(Some))
}

fn main() -> Result<()> {
    env_logger::init();
    let actions = Actions::from_args();
    match actions {
        Actions::GenerateParams {
            lower_limit,
            use_bigint,
        } => {
            let (prime, generator) = if use_bigint {
                let lower_limit: BigUint = pi(&lower_limit)?;
                let prime = generate_prime(lower_limit)?;
                debug!("Prime number is {}", prime);
                let generator = find_primitive_root(prime.clone());
                (prime.to_string(), generator.map(|g| g.to_string()))
            } else {
                let lower_limit = u128::from_str_radix(&lower_limit, 10)
                    .map_err(|_| Error::msg("Invalid lower_limit"))?;
                let prime = generate_prime(lower_limit).expect("result with u128 limit");
                debug!("Prime number is {}", prime);
                let generator = find_primitive_root(prime.clone());
                (prime.to_string(), generator.map(|g| g.to_string()))
            };
            match generator {
                Some(g) => {
                    println!("Prime number is {}", prime);
                    println!("And it's primitive root is {}", g);
                }
                None => {
                    eprintln!("ERROR: did not find prim. root");
                    exit(1)
                }
            }
        }
        Actions::Permutation {
            multiplier,
            modulo,
            seed,
            upper_limit,
            take,
            use_bigint,
        } => {
            if use_bigint {
                let multiplier = pi(multiplier)?;
                let modulo: BigUint = pi(modulo)?;

                let seed = flip(seed.map(|s| pi(s)))?.unwrap_or_else(|| {
                    let mut rng = rand::thread_rng();
                    rng.gen_range(BigUint::one()..modulo.clone())
                });
                let limit = flip(upper_limit.map(|s| pi(s)))?.unwrap_or_else(|| modulo.clone());
                let take = take.unwrap_or(usize::MAX);

                if multiplier >= modulo || seed < BigUint::one() {
                    bail!("Invalid multiplier, must be positive and less then modulo")
                }

                if seed >= modulo || seed < BigUint::one() {
                    bail!("Invalid seed, must be positive and less then modulo")
                }

                if limit > modulo || seed <= BigUint::one() {
                    bail!("Invalid limit, must be bigger then 1 and less then  or equal to modulo")
                }

                let gen = Generator::with_limit_and_state(modulo, multiplier, seed, limit);
                gen.take(take)
                    .enumerate()
                    .for_each(|(idx, n)| println!("{}:{}", idx + 1, n));
            } else {
                let multiplier: u128 = pi(multiplier)?;
                let modulo: u128 = pi(modulo)?;
                let seed = flip(seed.map(|s| pi(s)))?.unwrap_or_else(|| {
                    let mut rng = rand::thread_rng();
                    rng.gen_range(1..modulo)
                });
                let limit = flip(upper_limit.map(|s| pi(s)))?.unwrap_or_else(|| modulo.clone());
                let take = take.unwrap_or(usize::MAX);

                if multiplier >= modulo || seed < 1 {
                    bail!("Invalid multiplier, must be positive and less then modulo")
                }

                if seed >= modulo || seed < 1 {
                    bail!("Invalid seed, must be positive and less then modulo")
                }

                if limit > modulo || seed <= 1 {
                    bail!("Invalid limit, must be bigger then 1 and less then  or equal to modulo")
                }

                let gen = Generator::with_limit_and_state(modulo, multiplier, seed, limit);
                gen.take(take)
                    .enumerate()
                    .for_each(|(idx, n)| println!("{}:{}", idx + 1, n));
            }
        }
        Actions::RandomPermutation {
            upper_limit,
            seed,
            take,
        } => {
            let seed = seed.map(|s| sha2::Sha256::digest(s.as_bytes()));
            let p = random_permutation(upper_limit, seed);
            let sz = p.len();
            p.into_iter()
                .take(take.unwrap_or(sz))
                .enumerate()
                .for_each(|(idx, i)| println!("{}:{}", idx, i))
        }
    }
    Ok(())
}
