import aiohttp
import asyncio
import pandas as pd
import time
import json
import os
from collections import defaultdict

# Just the core testing functions needed
def parse_aa_string(aa_string):
    """Parse strings like '2:H,V' into list of amino acids"""
    if not aa_string:
        return []
    count, aas = aa_string.split(':')
    return aas.split(',')

async def fetch_data(session, protein, position):
    """Fetch data for a single variant"""
    url = "https://alphamissense.hegelab.org/hotspotapi"
    params = {"uid": protein, "resi": position}
    try:
        async with session.get(url, params=params) as response:
            if response.status == 200:
                data = await response.json()
                # Use strict pathogenic classification
                data['pathogenic_aas'] = parse_aa_string(data.get('pathogenic', ''))
                return data
            else:
                print(f"Error: Status {response.status}")
    except Exception as e:
        print(f"Error during request: {e}")
    return None

async def test_variant():
    """Test with TMEM88B example"""
    # Create test data
    test_data = {
        'Protein': ['TMEM88B'],
        'Position': [19],
        'New_AA': ['N']
    }
    df = pd.DataFrame(test_data)
    
    print("Testing variant: TMEM88B D19N")
    
    async with aiohttp.ClientSession() as session:
        for _, row in df.iterrows():
            data = await fetch_data(session, row['Protein'], row['Position'])
            if data:
                print("\nAPI Response:")
                print(json.dumps(data, indent=2))
                
                print("\nClassification:")
                if row['New_AA'] in data['pathogenic_aas']:
                    print(f"Result: {row['New_AA']} is PATHOGENIC")
                else:
                    print(f"Result: {row['New_AA']} is NOT pathogenic")
                    print("Pathogenic amino acids:", data['pathogenic_aas'])

if __name__ == "__main__":
    asyncio.run(test_variant())