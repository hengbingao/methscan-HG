# GitHub Repository Setup Guide

## 📝 Step-by-Step Instructions

### 1. Create GitHub Repository

1. Go to https://github.com/new
2. Repository name: `methscan2`
3. Description: "Single-cell DNA methylation analysis toolkit"
4. Choose: Public or Private
5. **DO NOT** initialize with README, .gitignore, or license (we already have these)
6. Click "Create repository"

### 2. Initialize Local Git Repository

```bash
cd /path/to/methscan2-project

# Initialize git
git init

# Add all files
git add .

# First commit
git commit -m "Initial commit: MethSCAn2 v0.1.0"

# Add remote (replace 'yourusername' with your GitHub username)
git remote add origin https://github.com/yourusername/methscan2.git

# Push to GitHub
git branch -M main
git push -u origin main
```

### 3. Installation for Users

After pushing to GitHub, users can install with:

```bash
# Direct installation from GitHub
pip install git+https://github.com/yourusername/methscan2.git

# Or with specific branch/tag
pip install git+https://github.com/yourusername/methscan2.git@main
```

### 4. Development Installation

For contributors:

```bash
# Clone repository
git clone https://github.com/yourusername/methscan2.git
cd methscan2

# Install in editable mode
pip install -e .

# Install development dependencies
pip install -e ".[dev]"
```

### 5. Create a Release

When ready to release v0.1.0:

```bash
# Tag the release
git tag -a v0.1.0 -m "Release version 0.1.0"
git push origin v0.1.0
```

Then on GitHub:
1. Go to your repository
2. Click "Releases" → "Create a new release"
3. Select the v0.1.0 tag
4. Add release notes
5. Publish

### 6. (Optional) Publish to PyPI

To make installation even easier with `pip install methscan2`:

```bash
# Install build tools
pip install build twine

# Build distribution
python -m build

# Upload to PyPI (you'll need a PyPI account)
twine upload dist/*
```

Then users can simply:
```bash
pip install methscan2
```

## 🎯 Quick Test

After installation, test it works:

```python
import methscan2 as ms2
print(ms2.__version__)  # Should print: 0.1.0

# Create a simple test
import numpy as np
import pandas as pd

# Simulate data
X = np.random.rand(10, 100)
var = pd.DataFrame({
    'chrom': ['chr1'] * 100,
    'start': np.arange(100) * 100,
    'end': (np.arange(100) + 1) * 100
})

mdata = ms2.MethylationData(X=X, var=var)
print(f"Created: {mdata.n_obs} cells × {mdata.n_vars} sites")
```

## 📚 Documentation Options

### Option 1: GitHub Pages (Simple)

1. Create `docs/` folder with markdown files
2. Go to Settings → Pages
3. Source: Deploy from main branch, /docs folder

### Option 2: Read the Docs (Advanced)

1. Sign up at https://readthedocs.org
2. Import your GitHub repository
3. Add `docs/conf.py` for Sphinx configuration

## 🔧 Maintenance

### Update version

Edit `methscan2/__init__.py`:
```python
__version__ = "0.2.0"
```

And `setup.py`:
```python
version="0.2.0",
```

### Add new features

```bash
# Create feature branch
git checkout -b feature/new-analysis

# Make changes
# ... code ...

# Commit and push
git add .
git commit -m "Add new analysis feature"
git push origin feature/new-analysis

# Create Pull Request on GitHub
```

## ✅ Checklist

Before going public:

- [ ] All core functions implemented
- [ ] Tests passing (run `pytest`)
- [ ] README is clear and informative
- [ ] Example scripts work
- [ ] Version number set correctly
- [ ] License file included
- [ ] .gitignore configured
- [ ] GitHub repository created
- [ ] Code pushed to GitHub
- [ ] Installation tested from GitHub

## 🚀 You're Ready!

Your package is now:
- ✅ On GitHub
- ✅ Installable with pip
- ✅ Ready for collaboration
- ✅ Set up for CI/CD testing

Share your repository URL and users can start using MethSCAn2! 🎉
