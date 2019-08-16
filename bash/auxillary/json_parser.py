#!/Users/njwheeler/software/miniconda3/bin/python3

import asyncio
import requests
from concurrent.futures import ThreadPoolExecutor
# import json

# script, filename = argv


# WormBase ParaSite REST API documentation: https://parasite.wormbase.org/rest-13/documentation/info/homology_symbol
def fetch(session, id):
    paralogues = []
    print(id)

    # define the server
    server = "https://parasite.wormbase.org"
    # define the server extension and choose the condensed format with the paralogue type
    ext = "/rest-13/homology/symbol/caenorhabditis_elegans_PRJNA13758/" + \
        id + "?format=condensed;type=paralogues"
    # get the JSON response and decode the architecture to get the "id" entry; append to the list of paralogues and write it to files
    with session.get(server + ext, headers={"Content-Type": "application/json", "Accept": ""}) as response:
        decoded = response.json()
        for entry in decoded["data"][0]["homologies"]:
            paralogues.append(entry["id"])

    with open("/Users/njwheeler/GitHub/50HGI/auxillary/ChemoR/celegans_chemor_paralogues_geneid.txt", 'a') as f:
        f.write(id + "\n")
        paralogues = list(set(paralogues))

        for paralogue in paralogues:
            f.write(paralogue + "\n")


# for asynchronous server requests
# written using https://hackernoon.com/how-to-run-asynchronous-web-requests-in-parallel-with-python-3-5-without-aiohttp-264dc0f8546
# as a tutorial
async def get_data_async():

    filename = "/Users/njwheeler/GitHub/50HGI/auxillary/ChemoR/celegans_chemor_geneid.txt"
    gene_id = []

    with open(filename) as f:
        for line in f:
            line = line.strip(' \t\n\r')
            gene_id.append(line)

    with ThreadPoolExecutor(max_workers=10) as executor:
        with requests.Session() as session:

            # Initialize the event loop
            loop = asyncio.get_event_loop()

            # Use list comprehension to create a list of
            # tasks to complete. The executor will run the `fetch`
            # function for each csv in the csvs_to_fetch list
            tasks = [
                loop.run_in_executor(
                    executor,
                    fetch,
                    # Allows us to pass in multiple arguments to `fetch`
                    *(session, id)
                )
                for id in gene_id
            ]

            for response in await asyncio.gather(*tasks):
                pass


def main():
    loop = asyncio.get_event_loop()
    future = asyncio.ensure_future(get_data_async())
    loop.run_until_complete(future)

    paralogue_list = []

    with open("/Users/njwheeler/GitHub/50HGI/auxillary/ChemoR/celegans_chemor_paralogues_geneid.txt", 'r') as f:
        for line in f:
            line = line = line.strip(' \t\n\r')
            paralogue_list.append(line)

    paralogue_list = list(set(paralogue_list))
    with open("/Users/njwheeler/GitHub/50HGI/auxillary/ChemoR/celegans_chemor_paralogues_geneid2.txt", 'w') as f:

        for id in paralogue_list:
            f.write(id + "\n")


main()
