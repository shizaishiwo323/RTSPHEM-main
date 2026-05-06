#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GIF/MP4 Generation Helper Script
Called by MATLAB automation framework
"""

import os
import sys
import re
import argparse

# Fix Windows console encoding issue
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

from PIL import Image

def handle_long_path(path):
    """Handle Windows long path (>260 chars) by adding extended path prefix"""
    if sys.platform == 'win32' and len(path) > 200:
        # Convert to absolute path and add extended path prefix
        abs_path = os.path.abspath(path)
        if not abs_path.startswith('\\\\?\\'):
            return '\\\\?\\' + abs_path
    return path

def get_inversion_steps(inversion_folder):
    """Get valid step numbers from inversion_results folder"""
    steps = set()
    if not os.path.exists(inversion_folder):
        return steps
    
    for filename in os.listdir(inversion_folder):
        match = re.search(r'_t(\d{4})_', filename)
        if match:
            steps.add(int(match.group(1)))
    return sorted(steps)

def create_gif(input_folder, output_file, speed=1.0, resize=None, valid_steps=None):
    """Create GIF from images in folder"""
    images = []
    input_folder = handle_long_path(input_folder)
    filenames = sorted([f for f in os.listdir(input_folder) 
                       if f.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif'))])
    
    for filename in filenames:
        if 'timestep' in filename:
            match = re.search(r'timestep_(\d{4})', filename)
            if match:
                step_num = int(match.group(1))
                if valid_steps is not None and step_num not in valid_steps:
                    continue
            
            filepath = handle_long_path(os.path.join(input_folder, filename))
            img = Image.open(filepath)
            if resize:
                img = img.resize(resize)
            images.append(img)
    
    if not images:
        print(f"Warning: No matching images found in {input_folder}")
        return False
    
    duration = int(200 / speed)
    output_file = handle_long_path(output_file)
    images[0].save(
        output_file,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=0
    )
    print(f"GIF created: {output_file}, {len(images)} frames")
    return True

def create_gif_from_inversion(inversion_folder, output_file, speed=1.0, resize=None):
    """Create GIF from inversion_results folder images"""
    images = []
    inversion_folder = handle_long_path(inversion_folder)
    filenames = sorted([f for f in os.listdir(inversion_folder) 
                       if f.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.gif'))])
    
    for filename in filenames:
        filepath = handle_long_path(os.path.join(inversion_folder, filename))
        img = Image.open(filepath)
        if resize:
            img = img.resize(resize)
        images.append(img)
    
    if not images:
        print(f"Warning: No images found in {inversion_folder}")
        return False
    
    duration = int(200 / speed)
    output_file = handle_long_path(output_file)
    images[0].save(
        output_file,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=0
    )
    print(f"GIF created: {output_file}, {len(images)} frames")
    return True

def process_dissolution_folder(base_folder, speed=1.0, resize=None, output_format='gif', valid_steps_str=None):
    """Process dissolution results folder and generate GIFs"""
    base_folder = handle_long_path(base_folder)
    if not os.path.exists(base_folder):
        print(f"Error: Folder does not exist: {base_folder}")
        return False
    
    inversion_folder = handle_long_path(os.path.join(base_folder, 'inversion_results'))
    
    # Parse valid steps
    if valid_steps_str:
        valid_steps = set(int(x) for x in valid_steps_str.split(',') if x.strip())
    else:
        valid_steps = get_inversion_steps(inversion_folder)
    
    if not valid_steps:
        print("Warning: No valid steps found")
        valid_steps = None
    else:
        print(f"Valid steps: {sorted(valid_steps)}")
    
    ext = output_format.lower()
    success = True
    
    # Generate timestep GIF
    timestep_output = handle_long_path(os.path.join(base_folder, f"animation_timestep.{ext}"))
    if not create_gif(base_folder, timestep_output, speed, resize, valid_steps):
        success = False
    
    # Generate inversion GIF
    if os.path.exists(inversion_folder):
        inversion_output = handle_long_path(os.path.join(base_folder, f"animation_inversion.{ext}"))
        if not create_gif_from_inversion(inversion_folder, inversion_output, speed, resize):
            success = False
    
    return success

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate GIF animation')
    parser.add_argument('folder', help='Dissolution results folder path')
    parser.add_argument('--speed', type=float, default=0.5, help='Playback speed')
    parser.add_argument('--format', default='gif', help='Output format (gif/mp4)')
    parser.add_argument('--valid-steps', default=None, help='Valid step list, comma separated')
    
    args = parser.parse_args()
    
    success = process_dissolution_folder(
        args.folder,
        speed=args.speed,
        output_format=args.format,
        valid_steps_str=getattr(args, 'valid_steps', None)
    )
    
    sys.exit(0 if success else 1)
