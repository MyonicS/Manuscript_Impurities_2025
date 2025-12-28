
import pandas as pd
import numpy as np
import pathlib as Path

def parse_xrd_file(file_path):
    """
    Parse Bruker XRD .uxd file and extract metadata and measurement data.
    
    Parameters:
    -----------
    file_path : str
        Path to the .uxd file
    
    Returns:
    --------
    dict : Dictionary containing:
        - 'data': pandas DataFrame with '2theta' and 'counts' columns
        - 'metadata': dict with measurement parameters
        - 'sample_name': str with sample identifier
    """
    
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
            lines = file.readlines()
        
        # Initialize containers
        metadata = {}
        data_lines = []
        data_started = False
        
        # Parse file line by line
        for line in lines:
            line = line.strip()
            
            # Check for data header (indicates start of measurement data)
            if '2THETA' in line and 'Cnt2_D1' in line and line.startswith(';'):
                data_started = True
                continue
                
            # Skip empty lines and other comments
            if not line or (line.startswith(';') and not data_started):
                continue
            
            # If we're in the data section, collect the measurement points
            if data_started:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        # Parse 2theta and counts values
                        two_theta = float(parts[0])
                        counts = float(parts[1])
                        data_lines.append([two_theta, counts])
                    except (ValueError, IndexError):
                        # Skip lines that can't be parsed as data
                        continue
            
            # Parse metadata parameters (lines starting with _)
            elif line.startswith('_'):
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip().replace('_', '')
                    value = value.strip()
                    metadata[key] = value
        
        # Create DataFrame from measurement data
        if data_lines:
            df = pd.DataFrame(data_lines, columns=['2theta', 'counts'])
            df = df.sort_values('2theta').reset_index(drop=True)
            
            # Remove any rows with NaN or infinite values
            df = df.replace([np.inf, -np.inf], np.nan).dropna()
            
            # Check if we still have data after cleaning
            if len(df) == 0:
                print(f"⚠️ Warning: No valid measurement data after cleaning in {file_path}")
                return None
        else:
            print(f"⚠️ Warning: No measurement data found in {file_path}")
            return None
        
        # Extract sample name from metadata or filename
        sample_name = metadata.get('SAMPLE', '')
        if not sample_name:
            sample_name = metadata.get('V4_INF_SAMPLEID', '')
        if not sample_name:
            sample_name = Path(file_path).stem
        
        result = {
            'data': df,
            'metadata': metadata,
            'sample_name': sample_name,
            'file_path': file_path
        }
        
        return result
    
    except Exception as e:
        print(f"❌ Error parsing {file_path}: {str(e)}")
        return None

