#!/usr/bin/env python3
"""
Quick ALLC to COV Conversion

Minimal script that replicates your bash workflow:
1. Convert ALLC → COV
2. Separate CG, CH contexts
3. Calculate global methylation

Usage:
    python convert_allc.py --input raw/allc --output processed --jobs 100
"""

import argparse
import methscan2 as ms2


def main():
    parser = argparse.ArgumentParser(
        description='Convert ALLC files to COV format and calculate global methylation'
    )
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input directory containing .allc.tsv.gz files'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='processed',
        help='Output directory (default: processed)'
    )
    
    parser.add_argument(
        '--contexts', '-c',
        nargs='+',
        default=['CG', 'CH', 'all'],
        help='Contexts to extract (default: CG CH all)'
    )
    
    parser.add_argument(
        '--jobs', '-j',
        type=int,
        default=-1,
        help='Number of parallel jobs (default: -1 = all CPUs)'
    )
    
    parser.add_argument(
        '--no-global',
        action='store_true',
        help='Skip global methylation calculation'
    )
    
    args = parser.parse_args()
    
    # Run pipeline
    print(f"\nProcessing ALLC files from: {args.input}")
    print(f"Output directory: {args.output}")
    print(f"Contexts: {args.contexts}")
    print(f"Parallel jobs: {args.jobs if args.jobs > 0 else 'all CPUs'}")
    
    global_stats = ms2.pp.batch_process_allc(
        allc_dir=args.input,
        output_base_dir=args.output,
        contexts=args.contexts,
        calculate_global=not args.no_global,
        n_jobs=args.jobs,
        verbose=True
    )
    
    if not args.no_global:
        print("\n" + "=" * 60)
        print("Global Methylation Summary")
        print("=" * 60)
        
        for ctx, df in global_stats.items():
            mean_level = df[f'global_{ctx}_level'].mean()
            mean_sites = df[f'{ctx}_sites'].mean()
            print(f"\n{ctx}:")
            print(f"  Mean methylation: {mean_level:.3f}")
            print(f"  Mean sites/cell:  {mean_sites:.0f}")
    
    print("\n✓ Complete!")


if __name__ == '__main__':
    main()
