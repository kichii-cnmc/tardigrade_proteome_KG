#!/bin/bash

# Shared utility for logging download metadata
# Source this file in your download scripts with: source "$(dirname "$0")/download_logger.sh"

log_download() {
    local file="$1"
    local url="$2"
    local source="$3"
    local organism="$4"
    local taxonomy_id="$5"
    local version="$6"
    local script_name="$(basename "$0")"
    local timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    
    # Get project root directory (assuming we're in src/1_collector/)
    local project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    mkdir -p "$project_root/logs"
    local log_file="$project_root/logs/download_metadata.json"
    
    echo "Logging download: $organism - $source"
    
    # Initialize file if it doesn't exist
    if [ ! -f "$log_file" ]; then
        echo '{"downloads": []}' > "$log_file"
    fi
    
    # Create temporary file for JSON manipulation
    local temp_file=$(mktemp)
    
    # Create new entry (escape quotes in JSON values)
    local escaped_url=$(echo "$url" | sed 's/"/\\"/g')
    local escaped_source=$(echo "$source" | sed 's/"/\\"/g')
    local escaped_organism=$(echo "$organism" | sed 's/"/\\"/g')
    local escaped_file=$(echo "$file" | sed 's/"/\\"/g')
    
    cat << EOF > "$temp_file"
{
  "timestamp": "$timestamp",
  "script": "$script_name",
  "source": "$escaped_source",
  "file": "$escaped_file",
  "url": "$escaped_url",
  "organism": "$escaped_organism",
  "taxonomy_id": "$taxonomy_id",
  "version": "$version"
}
EOF
    
    # Add entry to JSON array using Python (more reliable than sed)
    if command -v python3 &> /dev/null; then
        python3 -c "
import json
import sys

# Read existing log file
try:
    with open('$log_file', 'r') as f:
        data = json.load(f)
except:
    data = {'downloads': []}

# Read new entry
with open('$temp_file', 'r') as f:
    new_entry = json.load(f)

# Add new entry
data['downloads'].append(new_entry)

# Write back
with open('$log_file', 'w') as f:
    json.dump(data, f, indent=2)
"
    else
        # Fallback: simple append (less reliable but works without Python)
        if grep -q '"downloads": \[\]' "$log_file"; then
            # Empty array - replace with first entry
            sed -i '' "s/\"downloads\": \[\]/\"downloads\": [$(cat "$temp_file")]/" "$log_file"
        else
            # Has entries - append to array (this is fragile but better than nothing)
            echo "Warning: Using basic JSON append - install python3 for reliable JSON handling"
            sed -i '' "s/\]\s*$/,$(cat "$temp_file")]/" "$log_file"
        fi
    fi
    
    rm -f "$temp_file"
    echo "Download logged to $log_file"
}

# Function to show recent downloads
show_recent_downloads() {
    local project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    local log_file="$project_root/logs/download_metadata.json"
    
    if [ -f "$log_file" ]; then
        echo "Recent downloads:"
        # Simple display of the last few downloads (requires jq for proper parsing)
        if command -v jq &> /dev/null; then
            jq -r '.downloads[-5:] | .[] | "\(.timestamp) - \(.organism) from \(.source)"' "$log_file"
        else
            echo "Install 'jq' for formatted output, or check $log_file directly"
        fi
    else
        echo "No download history found"
    fi
}