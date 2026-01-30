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

# Add new entry at the beginning (newest first)
data['downloads'].insert(0, new_entry)

# Write back
with open('$log_file', 'w') as f:
    json.dump(data, f, indent=2)
"
    else
        # Fallback: simple prepend (less reliable but works without Python)
        if grep -q '"downloads": \[\]' "$log_file"; then
            # Empty array - replace with first entry
            sed -i '' "s/\"downloads\": \[\]/\"downloads\": [$(cat "$temp_file")]/" "$log_file"
        else
            # Has entries - prepend to array
            echo "Warning: Using basic JSON prepend - install python3 for reliable JSON handling"
            sed -i '' "s/\"downloads\": \[/\"downloads\": [$(cat "$temp_file"),/" "$log_file"
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
        echo "Recent downloads (newest first):"
        # Simple display of the first few downloads (requires jq for proper parsing)
        if command -v jq &> /dev/null; then
            jq -r '.downloads[0:5] | .[] | "\(.timestamp) - \(.organism) from \(.source)"' "$log_file"
        else
            echo "Install 'jq' for formatted output, or check $log_file directly"
        fi
    else
        echo "No download history found"
    fi
}

# Function to archive old logs
archive_logs() {
    local project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    local log_file="$project_root/logs/download_metadata.json"
    local archive_dir="$project_root/logs/archive"
    
    if [ -f "$log_file" ]; then
        mkdir -p "$archive_dir"
        local archive_name="download_metadata_$(date +%Y%m%d_%H%M%S).json"
        mv "$log_file" "$archive_dir/$archive_name"
        echo "Logs archived to $archive_dir/$archive_name"
        echo '{"downloads": []}' > "$log_file"
        echo "Fresh log file created"
    else
        echo "No log file to archive"
    fi
}

# Function to clean old logs (keep only last N entries)
clean_logs() {
    local keep_count=${1:-50}  # Default: keep last 50 entries
    local project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    local log_file="$project_root/logs/download_metadata.json"
    
    if [ -f "$log_file" ] && command -v python3 &> /dev/null; then
        python3 -c "
import json

with open('$log_file', 'r') as f:
    data = json.load(f)

# Keep only the first $keep_count entries (newest)
data['downloads'] = data['downloads'][:$keep_count]

with open('$log_file', 'w') as f:
    json.dump(data, f, indent=2)

print('Kept most recent $keep_count download entries')
"
    else
        echo "Python3 required for log cleaning"
    fi
}