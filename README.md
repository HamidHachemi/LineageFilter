# LineageFilter

LineageFilter is a Python library enabling the estimation of the taxonomic composition of complex samples using metaproteomics and machine learning.

## Install
After downloading the package file (.whl file), install it as follows:
```bash
#Unix/macOS
python3 -m pip install --force-reinstall .\lineagefilter-0.0.1-py3-none-any.whl
#Windows
py -m pip install --force-reinstall .\lineagefilter-0.0.1-py3-none-any.whl
```

## Simple Demo

```python
from github import Github

# Authentication is defined via github.Auth
from github import Auth

# using an access token
auth = Auth.Token("access_token")

# First create a Github instance:

# Public Web Github
g = Github(auth=auth)

# Github Enterprise with custom hostname
g = Github(base_url="https://{hostname}/api/v3", auth=auth)

# Then play with your Github objects:
for repo in g.get_user().get_repos():
    print(repo.name)

# To close connections after use
g.close()
```
