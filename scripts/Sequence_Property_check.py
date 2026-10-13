import peptides

def calculate_sequence_metrics(sequence, angle):
    """
    Calculates the exact biochemical parameters used in S3-Dock Phase 1.
    """
    pep = peptides.Peptide(sequence)
    c = pep.charge(pH=7.4)
    b = pep.boman()
    m = pep.hydrophobic_moment(angle=angle)
    return c, b, m

def main():
    print("\n====================================================")
    print("🔬 S3-DOCK INTERACTIVE METRICS CHECKER")
    print("====================================================")
    print("Type 'exit' or 'quit' at any prompt to close the tool.\n")

    while True:
        # 1. Ask for the sequence
        seq = input("Enter amino acid sequence: ").strip().upper()
        
        if seq.lower() in ['exit', 'quit']:
            break
            
        if not seq:
            print("❌ Sequence cannot be empty. Try again.\n")
            continue

        # 2. Ask for the structure with a fast default
        struct_input = input("Target structure [helix/beta] (Press Enter for helix): ").strip().lower()
        
        if struct_input in ['exit', 'quit']:
            break
            
        if struct_input in ["beta", "beta_sheet", "b"]:
            calc_angle = 180
            struct_name = "Beta-Sheet"
        else:
            calc_angle = 100
            struct_name = "Alpha-Helix"

        # 3. Calculate and display
        try:
            charge, boman, moment = calculate_sequence_metrics(seq, calc_angle)
            
            print("\n----------------------------------------------------")
            print(f"Sequence:            {seq}")
            print(f"Assumed Structure:   {struct_name} (Angle: {calc_angle}°)")
            print(f"Net Charge (pH 7.4): {charge:.3f}")
            print(f"Boman Index:         {boman:.3f}")
            print(f"Hydrophobic Moment:  {moment:.3f}")
            print("----------------------------------------------------\n")
            
        except ValueError as e:
            print(f"\n❌ [ERROR] Invalid amino acid sequence. Make sure you are using standard 1-letter codes.\n")

if __name__ == "__main__":
    main()