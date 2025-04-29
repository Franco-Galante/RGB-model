# extracts the datasets from the 7z file(s) stored in the 'data' folder

import os
import py7zr


# extract the xml file from the 7z archives
def extract_7z(archive_path, destination_path):
    with py7zr.SevenZipFile(archive_path, mode='r') as z:
        z.extractall(path=destination_path)


if __name__ == '__main__':

    dataset         = str(os.environ.get('DATASET', 'math')) # 'cs' or 'poker'
    raw_data_folder = str(os.environ.get('RAW_DATA_FOLDER', 'math.stackexchange.com'))
    
    file_format = str(os.environ.get('OUT_FORMAT', 'csv')) # 'csv' or 'pkl'

    width = 80 # terminal width

    print('\n', '*** START EXTRACTING FILES ***'.center(width, '*'), '\n', flush=True)
    
    data_folder = os.path.join('data', raw_data_folder)

    out_posts = os.path.join(data_folder, 'Posts.xml')
    out_votes = os.path.join(data_folder, 'Votes.xml')
    out_badge = os.path.join(data_folder, 'Badges.xml')

    if not os.path.isdir(data_folder):
        os.makedirs(data_folder)

    if dataset == 'cs': # one 7z for each sub-dataset

        if not os.path.exists(out_posts):
            print('... extracting posts 7zip data ...'.center(width), flush=True)
            extract_7z(os.path.join('data','stackoverflow.com-Posts.7z'), data_folder)
            print('... xml posts file extracted ...'.center(width), flush=True)
        else:
            print('... xml posts file already extracted ...'.center(width), flush=True)
        
        if not os.path.exists(out_votes):
            print('... extracting votes 7zip data ...'.center(width), flush=True)
            extract_7z(os.path.join('data','stackoverflow.com-Votes.7z'), data_folder)
            print('... xml votes file extracted ...'.center(width), flush=True)
        else:
            print('... xml votes file already extracted ...'.center(width), flush=True)

        if not os.path.exists(out_badge):
            print('... extracting badge 7zip data ...'.center(width), flush=True)
            extract_7z(os.path.join('data','stackoverflow.com-Badges.7z'), data_folder)
            print('... xml badge file extracted ...'.center(width), flush=True)
        else:
            print('... xml badge file already extracted ...'.center(width), flush=True)

    else: # all sub-datasets in the same 7z

        if not os.path.exists(out_posts) or not os.path.exists(out_votes) or not os.path.exists(out_badge):
            print('... extracting 7zip data ...'.center(width), flush=True)
            extract_7z(os.path.join('data', raw_data_folder + '.7z'), data_folder)
            print('... xml  files extracted ...'.center(width), flush=True)
        else:
            print('... xml files already extracted ...'.center(width), flush=True)

    print('\n', '*** FILES EXTRACTED ***'.center(width, '*'), '\n', flush=True)
