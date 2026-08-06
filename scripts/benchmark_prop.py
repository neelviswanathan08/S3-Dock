import argparse
import peptides

def calculate_sequence_metrics(sequence, angle):
    """
    Calculates the exact biochemical parameters used in S3-Dock Phase 1.
    """
    # Initialize the peptide object
    pep = peptides.Peptide(sequence)
    
    # Calculate metrics matching your pipeline's logic
    c = pep.charge(pH=7.4)
    b = pep.boman()
    m = pep.hydrophobic_moment(angle=angle)
    
    return c, b, m

if __name__ == "__main__":
    # Set up command line arguments
    parser = argparse.ArgumentParser(description="Calculate sequence metrics (Charge, Boman, Hydrophobic Moment)")
    parser.add_argument("sequence", type=str, help="The amino acid sequence to evaluate (e.g., SQETASDAAKAL)")
    parser.add_argument("--structure", type=str, choices=["helix", "beta_sheet", "other"], default="helix",
                        help="Target secondary structure (determines angle for hydrophobic moment: helix=100, beta_sheet=180)")
    
    args = parser.parse_args()
    seq = args.sequence.upper()
    
    # Determine the calculation angle based on your pipeline's logic
    if args.structure == "helix":
        calc_angle = 100
    elif args.structure == "beta_sheet":
        calc_angle = 180
    else:
        calc_angle = 100

    try:
        # Run calculation
        charge, boman, moment = calculate_sequence_metrics(seq, calc_angle)
        
        # Print formatted output
        print("\n====================================================")
        print(f"🔬 S3-DOCK METRICS CHECKER")
        print("====================================================")
        print(f"Sequence:           {seq}")
        print(f"Assumed Structure:  {args.structure.capitalize()} (Angle: {calc_angle}°)")
        print("----------------------------------------------------")
        print(f"Net Charge (pH 7.4): {charge:.3f}")
        print(f"Boman Index:         {boman:.3f}")
        print(f"Hydrophobic Moment:  {moment:.3f}")
        print("====================================================\n")
        
    except ValueError as e:
        print(f"\n❌ [ERROR] Invalid amino acid sequence: {e}\n")