"""Example: Query a single star by name or alias.

Demonstrates how to use CorgiQuery to look up a star and inspect
basic stellar parameters returned by the corgidb API.
"""

from corgidb import CorgiQuery

cq = CorgiQuery()

# Query by primary name
df = cq.query_star("47 UMa")
print("Query: '47 UMa'")
print(df.to_string(index=False))
print()

# Another star
df2 = cq.query_star("eps Eri")
print("Query: 'eps Eri'")
print(df2.to_string(index=False))
print()

# Check for a star that does not exist — returns an empty DataFrame
df_missing = cq.query_star("NotAStar 9999")
print(f"Query: 'NotAStar 9999' — rows returned: {len(df_missing)}")
