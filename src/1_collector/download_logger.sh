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
    
    # Create new entry with proper JSON escaping
    local new_entry="{
      \"timestamp\": \"$timestamp\",
      \"script\": \"$script_name\",
      \"source\": \"$source\",
      \"file\": \"$file\",
      \"url\": \"$url\",
      \"organism\": \"$organism\",
      \"taxonomy_id\": \"$taxonomy_id\",
      \"version\": \"$version\"
    }"
    
    # Initialize file if it doesn't exist
    if [ ! -f "$log_file" ]; then
        echo '{"downloads": []}' > "$log_file"
    fi
    
    # Add new entry
    if grep -q '"downloads": \[\]' "$log_file"; then
        # Empty array - add first entry
        sed -i '' 's/"downloads": \[\]/"downloads": [\
    '"$new_entry"'\
  ]/' "$log_file"
    else
        # Has entries - add to array
        sed -i '' 's/\(.*\)\]/\1,\
    '"$new_entry"'\
  ]/' "$log_file"
    fi
    
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