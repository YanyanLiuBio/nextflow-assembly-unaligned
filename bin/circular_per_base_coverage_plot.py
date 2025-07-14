#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches

# Updated color scheme to match the target image
FEATURE_COLORS = {
    'gene': '#4472C4', 'CDS': '#70AD47', 'promoter': '#FFC000',
    'terminator': '#C5504B', 'origin': '#7030A0', 'rep_origin': '#7030A0',
    'misc_feature': '#7F7F7F', 'polyA_signal': '#8B4513',
    'protein_bind': '#00B0F0', 'enhancer': '#FF69B4'
}



def create_circular_plasmid_plot(genbank_file, coverage_file, output_file, show_variants=False):
    gb_record = SeqIO.read(genbank_file, "genbank")
    plasmid_length = len(gb_record.seq)

    full_coverage_data = pd.read_csv(coverage_file)
    coverage_data = pd.DataFrame({
        'pos': full_coverage_data['Position'],
        'coverage': full_coverage_data['TotalReads']
    })

    high_mismatch_positions = []
    if show_variants and 'Mismatches' in full_coverage_data.columns:
        mismatch_rate = full_coverage_data['Mismatches'] / full_coverage_data['TotalReads'].replace(0, np.nan)
        high_mismatch_positions = full_coverage_data['Position'][mismatch_rate > 0.1].values

    # Create figure with white background
    fig = plt.figure(figsize=(14, 12), facecolor='white')  # Made wider to accommodate legend
    
    # Create a grid layout: main plot on left, legend on right
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.1)
    ax = fig.add_subplot(gs[0], projection='polar', facecolor='white')
    legend_ax = fig.add_subplot(gs[1], facecolor='white')
    
    # Set up angles for coverage data
    angles = 2 * np.pi * (coverage_data['pos'] - 1) / plasmid_length
    
    # Plot coverage with light blue fill, more transparent
    ax.fill_between(angles, 0, coverage_data['coverage'], alpha=0.4, color='lightblue', linewidth=0)
    
    # Add coverage outline for better definition
    ax.plot(angles, coverage_data['coverage'], color='lightblue', alpha=0.6, linewidth=0.5)

    # Highlight high mismatch positions if requested
    for pos in high_mismatch_positions:
        angle = 2 * np.pi * (pos - 1) / plasmid_length
        cov = coverage_data[coverage_data['pos'] == pos]['coverage'].values
        if len(cov) > 0:
            ax.plot([angle, angle], [0, cov[0]], color='red', linewidth=2, alpha=0.7)

    max_cov = coverage_data['coverage'].max()
    
    # Position features outside the coverage area
    feature_radius = max_cov * 1.05
    
    # Draw the black circular outline
    circle_angles = np.linspace(0, 2*np.pi, 1000)
    ax.plot(circle_angles, [feature_radius] * 1000, color='black', linewidth=3)

    # Track which feature types are actually present for the legend
    present_features = set()

    # Process and draw features
    for feature in gb_record.features:
        if feature.type in FEATURE_COLORS:
            present_features.add(feature.type)
            start, end = int(feature.location.start), int(feature.location.end)
            
            # Handle circular features that wrap around
            if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
                # Handle compound locations (features that span the origin)
                for part in feature.location.parts:
                    part_start, part_end = int(part.start), int(part.end)
                    start_angle = 2 * np.pi * part_start / plasmid_length
                    end_angle = 2 * np.pi * part_end / plasmid_length
                    
                    color = FEATURE_COLORS.get(feature.type, 'gray')
                    arc = np.linspace(start_angle, end_angle, max(10, int((part_end - part_start) / 10)))
                    ax.plot(arc, [feature_radius] * len(arc), color=color, linewidth=12, alpha=0.9, solid_capstyle='butt')
            else:
                start_angle = 2 * np.pi * start / plasmid_length
                end_angle = 2 * np.pi * end / plasmid_length
                
                color = FEATURE_COLORS.get(feature.type, 'gray')
                
                # Handle features that cross the origin (0/2π boundary)
                if end < start:  # Feature wraps around
                    # Draw first part (start to end of plasmid)
                    arc1 = np.linspace(start_angle, 2*np.pi, max(10, int((plasmid_length - start) / 10)))
                    ax.plot(arc1, [feature_radius] * len(arc1), color=color, linewidth=12, alpha=0.9, solid_capstyle='butt')
                    
                    # Draw second part (beginning of plasmid to end)
                    arc2 = np.linspace(0, end_angle, max(10, int(end / 10)))
                    ax.plot(arc2, [feature_radius] * len(arc2), color=color, linewidth=12, alpha=0.9, solid_capstyle='butt')
                else:
                    arc = np.linspace(start_angle, end_angle, max(10, int((end - start) / 10)))
                    ax.plot(arc, [feature_radius] * len(arc), color=color, linewidth=12, alpha=0.9, solid_capstyle='butt')

    # Add feature labels with leader lines - improved overlap prevention
    label_radius = feature_radius * 1.25
    all_features = []
    
    # First pass: collect all features with their positions and sizes
    for feature in gb_record.features:
        if feature.type in FEATURE_COLORS:
            start, end = int(feature.location.start), int(feature.location.end)
            
            # Calculate feature size first to filter out small features
            feature_size = end - start if end > start else (plasmid_length - start + end)
            
            # Skip labeling for features smaller than 50 bp
            if feature_size < 50:
                continue
            
            # Check match_length qualifier and skip if < 20%
            if 'match_length' in feature.qualifiers:
                try:
                    match_length = float(feature.qualifiers['match_length'][0])
                    if match_length < 20.0:
                        continue
                except (ValueError, IndexError):
                    # If match_length can't be parsed, proceed with labeling
                    pass
            
            # Get feature label
            label = None
            for key in ["label", "gene", "product", "note"]:
                if key in feature.qualifiers and feature.qualifiers[key][0].strip():
                    label = feature.qualifiers[key][0].strip()
                    break
            
            if not label:
                label = f"{feature.type}_{start}"
            
            if end > start:
                mid_pos = (start + end) / 2
            else:  # Wraps around origin
                mid_pos = ((start + end + plasmid_length) / 2) % plasmid_length
            
            mid_angle = 2 * np.pi * mid_pos / plasmid_length
            
            all_features.append({
                'label': label,
                'mid_angle': mid_angle,
                'size': feature_size,
                'start': start,
                'end': end
            })
    
    # Sort features by size (largest first) to prioritize important features
    all_features.sort(key=lambda x: x['size'], reverse=True)
    
    # Second pass: place labels with improved overlap avoidance
    placed_labels = []  # Store angle and label width for each placed label
    
    for feature_info in all_features:
        label = feature_info['label']
        mid_angle = feature_info['mid_angle']
        feature_size = feature_info['size']
        
        # Clean up label text
        clean_label = label.replace('_', ' ')
        if len(clean_label) > 30:
            clean_label = clean_label[:27] + "..."
        
        # Estimate label angular width based on text length and font size
        # Rough approximation: each character takes about 0.015 radians at the label radius
        char_width = 0.015 if feature_size > plasmid_length * 0.02 else 0.012
        label_angular_width = len(clean_label) * char_width
        
        # Find best position for this label with proper spacing
        best_angle = mid_angle
        min_overlap = float('inf')
        
        # Generate test angles in a wider range for better distribution
        test_offsets = [0]  # Start with preferred position
        for i in range(1, 25):  # Try more positions
            offset = i * np.pi / 36  # 5-degree increments
            test_offsets.extend([offset, -offset])
        
        for offset in test_offsets:
            test_angle = (mid_angle + offset) % (2 * np.pi)
            
            # Calculate overlap with all existing labels
            total_overlap = 0
            has_conflict = False
            
            for placed_angle, placed_width in placed_labels:
                # Calculate angular distance (shortest path around circle)
                angular_dist = min(abs(test_angle - placed_angle), 
                                 2*np.pi - abs(test_angle - placed_angle))
                
                # Required separation is half of each label's width plus small buffer
                required_separation = (label_angular_width + placed_width) / 2 + 0.05
                
                if angular_dist < required_separation:
                    overlap = required_separation - angular_dist
                    total_overlap += overlap
                    has_conflict = True
            
            # Prefer positions with no conflicts, otherwise minimize overlap
            if not has_conflict:
                best_angle = test_angle
                break
            elif total_overlap < min_overlap:
                best_angle = test_angle
                min_overlap = total_overlap
        
        # Store this label's position and width for future overlap calculations
        placed_labels.append((best_angle, label_angular_width))
        
        # Choose labeling style and radius based on feature size and conflicts
        label_radius_multiplier = 1.0
        
        # If there are still conflicts after our best effort, try different radii
        if min_overlap > 0 and len(placed_labels) > 8:  # If crowded, use multiple levels
            # Alternate between inner and outer label rings
            label_radius_multiplier = 0.9 if len(placed_labels) % 2 == 0 else 1.1
        
        current_label_radius = label_radius * label_radius_multiplier
        
        if feature_size <= plasmid_length * 0.02:  # Small features
            # Compact labels close to plasmid
            ax.annotate(clean_label, 
                       xy=(best_angle, feature_radius * 1.12 * label_radius_multiplier),
                       ha='center', va='center', 
                       fontsize=7, fontweight='normal',
                       bbox=dict(boxstyle="round,pad=0.2", 
                               facecolor="white", 
                               edgecolor='gray',
                               linewidth=0.3,
                               alpha=0.9),
                       annotation_clip=False,
                       zorder=10)
        else:
            # Larger features with leader lines
            # Draw leader line from feature midpoint to label
            ax.plot([mid_angle, best_angle], 
                   [feature_radius * 1.02, current_label_radius * 0.88], 
                   color='black', linewidth=0.8, alpha=0.7, zorder=5)
            
            # Add label with box background
            ax.annotate(clean_label, 
                       xy=(best_angle, current_label_radius),
                       ha='center', va='center', 
                       fontsize=9, fontweight='normal',
                       bbox=dict(boxstyle="round,pad=0.3", 
                               facecolor="white", 
                               edgecolor='black',
                               linewidth=0.5,
                               alpha=0.95),
                       annotation_clip=False,
                       zorder=10)

    # Add scale markers (kb markers)
    scale_radius = feature_radius * 1.4
    kb_positions = np.arange(0, plasmid_length, 1000)  # Every 1kb
    
    for kb_pos in kb_positions:
        if kb_pos < plasmid_length:
            angle = 2 * np.pi * kb_pos / plasmid_length
            kb_label = f"{kb_pos//1000} kb" if kb_pos > 0 else "0"
            
            # Draw scale line
            ax.plot([angle, angle], [feature_radius * 1.02, scale_radius * 0.95], 
                   color='gray', linewidth=1, alpha=0.6)
            
            # Add kb label
            ax.text(angle, scale_radius, kb_label, 
                   rotation=np.degrees(angle) - 90 if np.pi/2 < angle < 3*np.pi/2 else np.degrees(angle),
                   ha='center', va='center', fontsize=8, color='gray')

    # Configure the plot appearance
    ax.set_theta_zero_location('N')  # Start at top (12 o'clock)
    ax.set_theta_direction(-1)       # Clockwise
    ax.set_ylim(0, scale_radius * 1.1)
    
    # Explicitly draw the inner grid lines with darker color
    # Radial grid lines (spokes)
    num_spokes = 12  # 12 spokes like a clock
    for i in range(num_spokes):
        spoke_angle = 2 * np.pi * i / num_spokes
        ax.plot([spoke_angle, spoke_angle], [0, max_cov], 
               color='gray', linewidth=0.6, alpha=0.6, zorder=1)
    
    # Concentric circles at specific coverage levels (25, 50, 75, 100)
    coverage_levels = [25, 50, 75, 100]
    for cov_level in coverage_levels:
        if cov_level <= max_cov:  # Only draw if within the coverage range
            ax.plot(circle_angles, [cov_level] * 1000, 
                   color='gray', linewidth=0.6, alpha=0.6, zorder=1)
            
            # Add coverage level labels at 0 degrees (top)
            ax.text(0, cov_level, str(cov_level), 
                   ha='center', va='center', fontsize=8, 
                   color='gray', fontweight='bold',
                   bbox=dict(boxstyle="round,pad=0.2", 
                           facecolor="white", 
                           edgecolor='none',
                           alpha=0.8),
                   zorder=2)
    
    # Remove the default grid and ticks
    ax.set_rticks([])  # Remove radial tick labels
    ax.set_thetagrids([])  # Remove theta labels
    ax.grid(False)  # Turn off default grid since we drew our own
    ax.spines['polar'].set_visible(False)
    
    # Add outer circle near the kb labels
    outer_circle_radius = scale_radius * 0.85
    ax.plot(circle_angles, [outer_circle_radius] * 1000, color='black', linewidth=2)

    # Add title with exact length in base pairs
    base_name = os.path.splitext(os.path.basename(genbank_file))[0]
    title_name = base_name if gb_record.id.strip() in ["", "."] else gb_record.id
    
    ax.set_title(f"{title_name}\n{plasmid_length:,} bp", fontsize=14, fontweight='bold', pad=20)

    # Create legend for feature colors
    legend_ax.axis('off')  # Turn off axis for legend subplot
    
    # Create legend patches only for features actually present in the plasmid
    legend_patches = []
    feature_names = {
        'gene': 'Gene',
        'CDS': 'Coding Sequence',
        'promoter': 'Promoter',
        'terminator': 'Terminator',
        'origin': 'Origin',
        'rep_origin': 'Replication Origin',
        'misc_feature': 'Misc Feature',
        'polyA_signal': 'PolyA Signal',
        'protein_bind': 'Protein Binding',
        'enhancer': 'Enhancer'
    }
    
    # Sort present features for consistent legend order
    sorted_present_features = sorted(present_features, key=lambda x: list(FEATURE_COLORS.keys()).index(x))
    
    for feature_type in sorted_present_features:
        color = FEATURE_COLORS[feature_type]
        name = feature_names.get(feature_type, feature_type.replace('_', ' ').title())
        legend_patches.append(mpatches.Patch(color=color, label=name))
    
    # Add coverage legend item
    legend_patches.append(mpatches.Patch(color='lightblue', alpha=0.6, label='Coverage'))
    
    # Add high mismatch legend item if variants are shown
    if show_variants and len(high_mismatch_positions) > 0:
        legend_patches.append(mpatches.Patch(color='red', label='High Mismatch (>10%)'))
    
    # Position legend in the right subplot
    legend = legend_ax.legend(handles=legend_patches, 
                             loc='center left',
                             bbox_to_anchor=(0, 0.5),
                             fontsize=10,
                             frameon=True,
                             fancybox=True,
                             shadow=True,
                             title='Feature Types',
                             title_fontsize=12)
    
    # Set title font weight separately (matplotlib version compatibility)
    legend.get_title().set_fontweight('bold')

    # Save with high quality
    plt.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"✅ Saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Plot a circular plasmid map.")
    parser.add_argument("--pair-id", required=True, help="Sample or pair ID")
    parser.add_argument("--genbank", required=True, help="GenBank file (.gbk)")
    parser.add_argument("--coverage", required=True, help="Per-base CSV file")
    parser.add_argument("--variants", action="store_true", help="Highlight mismatches")

    args = parser.parse_args()

    output_file = f"{args.pair_id}_circular_coverage_plot.png"

    create_circular_plasmid_plot(
        genbank_file=args.genbank,
        coverage_file=args.coverage,
        output_file=output_file,
        show_variants=args.variants
    )

if __name__ == "__main__":
    main()
