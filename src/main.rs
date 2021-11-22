use std::process::exit;

use anyhow::{Error, Result, bail};
use log::debug;
use num_bigint::{BigUint};
use num_traits::{One};
use random_permutation::*;
use structopt::StructOpt;

#[derive(StructOpt)]
enum Actions {
    #[structopt(help="Generate params for Lehman LCG - prime for modulo and primitive root for multiplier")]
    GenerateParams {
        #[structopt(required = true, help = "Prime will be bigger then")]
        lower_limit: String,
        #[structopt(short = "b", long, help = "Use bigint - no upper limit, but slower (app. 7x)")]
        use_bigint: bool,
    },
    #[structopt(help="Prints psedorandom permutation")]
    Permutation {
        #[structopt(required = true, help = "Multiplier - must be primitive root of modulo parameter")]
        multiplier: String,
        #[structopt(required = true, help = "Modulo - must be prime")]
        modulo: String,
        #[structopt(short, long, help = "Seed value, must be smaller then modulo")]
        seed: Option<String>,
        #[structopt(short, long, help = "Upper limit for permutation, must be smaller then modulo")]
        upper_limit: Option<String>,
        #[structopt(short, long, help = "Print only first x numbers")]
        take: Option<usize> 
    },

    #[structopt(help="Prints random using Fisher-Yates and CSPRG")]
    RandomPermutation {
        #[structopt(required=true, help = "Max element in permutation - permutating 1..=this")]
        upper_limit: usize,
        #[structopt(short, long, help = "Seed value, must be smaller then modulo")]
        seed: Option<String>,
        #[structopt(short, long, help = "Print only first x numbers")]
        take: Option<usize> 
    }
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
                let lower_limit = bi(&lower_limit)?;
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
        Actions::Permutation { multiplier, modulo, seed, upper_limit , take} => {
            let multiplier = bi(multiplier)?;
            let modulo = bi(modulo)?;
            let seed = flip(seed.map(|s| bi(s)))?.unwrap_or_else(BigUint::one);
            let limit = flip(upper_limit.map(|s| bi(s)))?.unwrap_or_else(|| modulo.clone());
            let take = take.unwrap_or(usize::MAX);
            
            if multiplier>=modulo || seed < BigUint::one() {
                bail!("Invalid multiplier, must be positive and less then modulo")
            }

            if seed>=modulo || seed < BigUint::one() {
                bail!("Invalid seed, must be positive and less then modulo")
            }

            if limit>modulo || seed <= BigUint::one() {
                bail!("Invalid limit, must be bigger then 1 and less then  or equal to modulo")
            }

            let gen = Generator::with_limit_and_state(modulo, multiplier, seed, limit);
            gen.take(take).enumerate().for_each(|(idx,n)| println!("{}:{}", idx+1,n));
        },
        Actions::RandomPermutation { upper_limit, seed, take } => {
            let p = random_permutation(upper_limit);
            let sz = p.len();
            p.into_iter().take(take.unwrap_or(sz))
                .enumerate().for_each(|(idx,i)| println!("{}:{}", idx, i))
        }

        
    }
    Ok(())
}
