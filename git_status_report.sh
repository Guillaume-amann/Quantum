#!/bin/bash
# =============================================================================
#  git_status_report.sh — Comprehensive Git repository status visualization
#
#  Shows:
#    * Current branch and status
#    * All local branches with commit counts vs origin
#    * Uncommitted changes (staged, unstaged, untracked)
#    * Recent commits on each branch
#    * Branch tracking status (ahead/behind)
#    * Diverged branches needing rebase/merge
#
#  Usage:
#    chmod +x git_status_report.sh
#    ./git_status_report.sh
# =============================================================================

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# Check if we're in a git repo
if [ ! -d .git ]; then
    echo -e "${RED}✗ Not a git repository${NC}"
    exit 1
fi

REPO_NAME=$(basename "$(git rev-parse --show-toplevel)")

echo ""
echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}  Git Status Report: ${REPO_NAME}${NC}"
echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
echo ""

# =============================================================================
# Section 1: Current branch and HEAD status
# =============================================================================
echo -e "${CYAN}┌─ CURRENT BRANCH${NC}"

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
CURRENT_HASH=$(git rev-parse --short HEAD)
CURRENT_COMMIT=$(git log -1 --format="%s" HEAD)

echo -e "│ Branch: ${MAGENTA}${CURRENT_BRANCH}${NC}"
echo -e "│ Commit: ${GREEN}${CURRENT_HASH}${NC} — ${CURRENT_COMMIT}"
echo ""

# Get commit date
COMMIT_DATE=$(git log -1 --format="%ai" HEAD | cut -d' ' -f1,2)
echo -e "│ Date: ${COMMIT_DATE}"
echo -e "└─"
echo ""

# =============================================================================
# Section 2: Uncommitted changes
# =============================================================================
echo -e "${CYAN}┌─ UNCOMMITTED CHANGES${NC}"

STAGED=$(git diff --cached --name-only | wc -l)
UNSTAGED=$(git diff --name-only | wc -l)
UNTRACKED=$(git ls-files --others --exclude-standard | wc -l)

if [ $STAGED -eq 0 ] && [ $UNSTAGED -eq 0 ] && [ $UNTRACKED -eq 0 ]; then
    echo -e "│ ${GREEN}✓ Working tree clean${NC}"
else
    if [ $STAGED -gt 0 ]; then
        echo -e "│ ${GREEN}Staged:${NC} ${STAGED} file(s)"
        git diff --cached --name-only | sed 's/^/│   ✓ /'
    fi
    
    if [ $UNSTAGED -gt 0 ]; then
        echo -e "│ ${YELLOW}Modified:${NC} ${UNSTAGED} file(s)"
        git diff --name-only | sed 's/^/│   ⚡ /'
    fi
    
    if [ $UNTRACKED -gt 0 ]; then
        echo -e "│ ${RED}Untracked:${NC} ${UNTRACKED} file(s)"
        git ls-files --others --exclude-standard | sed 's/^/│   ? /'
    fi
fi

echo -e "└─"
echo ""

# =============================================================================
# Section 3: Branch tracking status
# =============================================================================
echo -e "${CYAN}┌─ TRACKING STATUS (Current branch: ${MAGENTA}${CURRENT_BRANCH}${CYAN})${NC}"

UPSTREAM=$(git rev-parse --abbrev-ref @{u} 2>/dev/null || echo "none")

if [ "$UPSTREAM" == "none" ]; then
    echo -e "│ ${YELLOW}⚠ No upstream set${NC}"
    echo -e "│ To set: git branch --set-upstream-to=origin/${CURRENT_BRANCH}"
else
    echo -e "│ Upstream: ${GREEN}${UPSTREAM}${NC}"
    
    # Get ahead/behind count
    AHEAD=$(git rev-list --count @{u}..HEAD 2>/dev/null || echo 0)
    BEHIND=$(git rev-list --count HEAD..@{u} 2>/dev/null || echo 0)
    
    if [ "$AHEAD" -eq 0 ] && [ "$BEHIND" -eq 0 ]; then
        echo -e "│ ${GREEN}✓ In sync${NC} (no commits ahead or behind)"
    elif [ "$AHEAD" -gt 0 ] && [ "$BEHIND" -eq 0 ]; then
        echo -e "│ ${YELLOW}↑ Ahead by ${AHEAD} commit(s)${NC} (ready to push)"
    elif [ "$BEHIND" -gt 0 ] && [ "$AHEAD" -eq 0 ]; then
        echo -e "│ ${YELLOW}↓ Behind by ${BEHIND} commit(s)${NC} (need to pull)"
    else
        echo -e "│ ${RED}✗ Diverged${NC}: Ahead ${AHEAD}, Behind ${BEHIND} (need rebase/merge)"
    fi
fi

echo -e "└─"
echo ""

# =============================================================================
# Section 4: All branches overview
# =============================================================================
echo -e "${CYAN}┌─ ALL BRANCHES${NC}"

# Get all branches with tracking info
git for-each-ref --format='%(refname:short)' refs/heads/ | while read branch; do
    BRANCH_HASH=$(git rev-parse --short "$branch")
    BRANCH_COMMIT=$(git log -1 --format="%s" "$branch")
    
    # Check if it's the current branch
    if [ "$branch" = "$CURRENT_BRANCH" ]; then
        PREFIX="${GREEN}→${NC}"
    else
        PREFIX=" "
    fi
    
    # Get upstream and sync status
    BRANCH_UPSTREAM=$(git rev-parse --abbrev-ref "$branch"@{u} 2>/dev/null || echo "")
    
    if [ -z "$BRANCH_UPSTREAM" ]; then
        SYNC_STATUS="${YELLOW}[no upstream]${NC}"
    else
        BRANCH_AHEAD=$(git rev-list --count "$branch"@{u}.."$branch" 2>/dev/null || echo 0)
        BRANCH_BEHIND=$(git rev-list --count "$branch".."$branch"@{u} 2>/dev/null || echo 0)
        
        if [ "$BRANCH_AHEAD" -eq 0 ] && [ "$BRANCH_BEHIND" -eq 0 ]; then
            SYNC_STATUS="${GREEN}[✓ sync]${NC}"
        elif [ "$BRANCH_AHEAD" -gt 0 ] && [ "$BRANCH_BEHIND" -eq 0 ]; then
            SYNC_STATUS="${YELLOW}[↑ ${BRANCH_AHEAD}]${NC}"
        elif [ "$BRANCH_BEHIND" -gt 0 ] && [ "$BRANCH_AHEAD" -eq 0 ]; then
            SYNC_STATUS="${YELLOW}[↓ ${BRANCH_BEHIND}]${NC}"
        else
            SYNC_STATUS="${RED}[↔ ${BRANCH_AHEAD}↑/${BRANCH_BEHIND}↓]${NC}"
        fi
    fi
    
    echo -e "│ ${PREFIX} ${MAGENTA}${branch}${NC} @ ${GREEN}${BRANCH_HASH}${NC} ${SYNC_STATUS}"
    echo -e "│   └─ ${BRANCH_COMMIT}"
done

echo -e "└─"
echo ""

# =============================================================================
# Section 5: Recent commit history (current branch)
# =============================================================================
echo -e "${CYAN}┌─ RECENT COMMITS (${MAGENTA}${CURRENT_BRANCH}${CYAN})${NC}"

git log --oneline -8 HEAD | nl | sed 's/^/│ /'

echo -e "└─"
echo ""

# =============================================================================
# Section 6: Stash status
# =============================================================================
echo -e "${CYAN}┌─ STASH${NC}"

STASH_COUNT=$(git stash list | wc -l)

if [ $STASH_COUNT -eq 0 ]; then
    echo -e "│ ${GREEN}No stashed changes${NC}"
else
    echo -e "│ ${YELLOW}${STASH_COUNT} stash(es)${NC}"
    git stash list | sed 's/^/│   /'
fi

echo -e "└─"
echo ""

# =============================================================================
# Section 7: Summary and recommendations
# =============================================================================
echo -e "${CYAN}┌─ SUMMARY & RECOMMENDATIONS${NC}"

NEEDS_ACTION=false

# Check for uncommitted changes
if [ $STAGED -gt 0 ] || [ $UNSTAGED -gt 0 ] || [ $UNTRACKED -gt 0 ]; then
    echo -e "│ ${YELLOW}⚠ Uncommitted changes:${NC}"
    [ $STAGED -gt 0 ] && echo -e "│   • Stage remaining changes: ${CYAN}git add .${NC}"
    [ $UNSTAGED -gt 0 ] && echo -e "│   • Commit changes: ${CYAN}git commit -m '...'${NC}"
    [ $UNTRACKED -gt 0 ] && echo -e "│   • Add/ignore untracked files: ${CYAN}git add <file>${NC}"
    NEEDS_ACTION=true
fi

# Check for tracking issues
if [ "$UPSTREAM" != "none" ]; then
    AHEAD=$(git rev-list --count @{u}..HEAD 2>/dev/null || echo 0)
    BEHIND=$(git rev-list --count HEAD..@{u} 2>/dev/null || echo 0)
    
    if [ "$AHEAD" -gt 0 ]; then
        echo -e "│ ${YELLOW}⚠ Commits to push:${NC} ${CYAN}git push origin ${CURRENT_BRANCH}${NC}"
        NEEDS_ACTION=true
    fi
    
    if [ "$BEHIND" -gt 0 ]; then
        echo -e "│ ${YELLOW}⚠ Updates to pull:${NC} ${CYAN}git pull origin ${CURRENT_BRANCH}${NC}"
        NEEDS_ACTION=true
    fi
    
    if [ "$AHEAD" -gt 0 ] && [ "$BEHIND" -gt 0 ]; then
        echo -e "│ ${RED}✗ Branch diverged:${NC}"
        echo -e "│   Option 1 (rebase): ${CYAN}git pull --rebase origin ${CURRENT_BRANCH}${NC}"
        echo -e "│   Option 2 (merge): ${CYAN}git pull origin ${CURRENT_BRANCH}${NC}"
        NEEDS_ACTION=true
    fi
fi

if [ ! "$NEEDS_ACTION" = true ]; then
    echo -e "│ ${GREEN}✓ All systems go!${NC} Ready to work."
fi

echo -e "└─"
echo ""

# =============================================================================
# Legend
# =============================================================================
echo -e "${CYAN}LEGEND${NC}"
echo -e "  ${GREEN}→${NC}  = Current branch"
echo -e "  ${GREEN}✓${NC}  = In sync with upstream"
echo -e "  ${YELLOW}↑${NC}  = Commits ahead (ready to push)"
echo -e "  ${YELLOW}↓${NC}  = Commits behind (need to pull)"
echo -e "  ${RED}✗${NC}  = Diverged (conflict resolution needed)"
echo ""

echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
echo ""