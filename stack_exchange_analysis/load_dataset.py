# Converts the xml into a Pandas df, filters it (to reduce its dimension) and 
# save the result into a csv (preferred) or pkl file for processing. The filters
# aren't strict by choice, further filtering can be performed in the processing 
# phase, according to the needs.


import os
import warnings
import itertools
import pandas as pd
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup


# reads a large file line by line in chunks (of 'chunk_size');
# the 'yield' keyword allows to return and later resume the function from where
# the execution left
def read_in_chunks(file_object, chunk_size):
    while True:
        lines = list(itertools.islice(file_object, chunk_size))
        if not lines: # no more lines in the file
            break
        yield lines


# the 'body' of the questions is in HTML format, we convert it to text. Lastly, 
# we extract the relevant information removing options and what is not needed 
# (e.g. the "Possible Duplicate" section).
def markdown_to_text(html_body):

    # the parsing may give warinings, catch and print them
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        soup = BeautifulSoup(html_body, "html.parser") # soup object (manipulate and
                                                    # navigate HTML object)
        if w:
            print(f"Warning: {w[-1].message}")
            print(f"Problematic object: {html_body}")


    # remove specific tags:
    # - blockquote (citing other posts)
    # - a (hyperlinks)
    # NOTE: so far I left the "edits" but may be removed
    tags_to_remove = ['blockquote', 'a', ]
    for tag in soup.find_all(tags_to_remove):
        tag.decompose() # remove the tag and its contents

    # convert in plain text (structure is preserved by the '\n' between objects)
    plain_text = soup.get_text().strip()
    
    return plain_text


if __name__ == '__main__':

    dataset         = str(os.environ.get('DATASET', 'math')) # 'cs' or 'poker'
    raw_data_folder = str(os.environ.get('RAW_DATA_FOLDER', 'math.stackexchange.com'))

    file_format = str(os.environ.get('OUT_FORMAT', 'csv')) # 'csv' or 'pkl'

    linesxchunk = int(os.environ.get('CHUNK_SIZE', 10000)) # lines to read in a chunk
    th_views    = int(os.environ.get('MIN_VIEWCOUNT', 0))  # min views on a question
    th_num_answ = int(os.environ.get('MIN_ANSWERS', 0))    # min number of answers to a question 

    no_checks    = False # debug flag
    progress_out = True  # print loading progress to stdout
    width = 80           # terminal width

    print('\n', '*** START PROCESSING DATASET ***'.center(width, '*'), '\n', flush=True)
    
    data_folder = os.path.join('data', raw_data_folder)

    out_posts = os.path.join(data_folder, 'Posts.xml')
    out_votes = os.path.join(data_folder, 'Votes.xml')
    out_badge = os.path.join(data_folder, 'Badges.xml')


    # create Pandas (filtered) dataframe from the large xml file

    # Filters over the 'Posts.xml' dataset:
    # - 'PostTypeId' \in ['1', '2'] (only questions and answers)
    # - [Q] int('AnswerCount') >= 2     answers had to compete
    # - [Q] int('Score') > 0            avoid ill-posed QUESTIONS
    # - [Q] int('ViewCount) > 0         arbitrary 'interest' threshold (none)

    out_folder = os.path.join('data', file_format, dataset)
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)

    # ------------------------- questions and answers --------------------------

    out_q   = os.path.join(out_folder, 'df_questions.' + file_format)
    out_a   = os.path.join(out_folder, 'df_answers.' + file_format)
    out_ids = os.path.join(out_folder, 'filtered_ids.csv') # always csv

    if not os.path.exists(out_q) or not os.path.exists(out_a) or not os.path.exists(out_ids):
        
        questions_row_list = [] # list of dictionaries
        answers_row_list   = []

        Q_attr_to_keep = ['Id', 'PostTypeId', 'AcceptedAnswerId', 'CreationDate', 'OwnerUserId', 'Title', 'Body']
        A_attr_to_keep = ['Id', 'PostTypeId', 'AcceptedAnswerId', 'CreationDate', 'Score', 'ViewCount', 'OwnerUserId', 'Title', 'ParentId']
        
        print('... saving Q&A to file ...'.center(width), flush=True)

        num_lines = sum(1 for _ in open(out_posts, 'rb'))
        lines_count, prev_prog = 0, -1 # textual progression (for server)
        single_answer_questions = 0

        with open(out_posts, 'r', encoding='utf-8') as file:
            for _ in range(2): # (dataset-specific) skip first two lines
                next(file)

            for chunk in read_in_chunks(file, linesxchunk): # read chunks 
                for line in chunk:
                    if line.strip() != '</posts>': # signals end of file
                        root = ET.fromstring(line) # parse XML line

                        if root.attrib.get('PostTypeId') == '1':   # question
                            row_dict = {}
                            if not pd.isna(root.attrib.get('Score')) and int(root.attrib.get('Score')) > 0:
                                row_dict['Score'] = int(root.attrib.get('Score'))
                            else:
                                continue # skip this row (line)

                            if not pd.isna(root.attrib.get('ViewCount')) and int(root.attrib.get('ViewCount')) > th_views:
                                row_dict['ViewCount'] = int(root.attrib.get('ViewCount'))
                            else:
                                continue

                            if not pd.isna(root.attrib.get('AnswerCount')) and int(root.attrib.get('AnswerCount')) >= th_num_answ:
                                row_dict['AnswerCount'] = int(root.attrib.get('AnswerCount'))
                            else:
                                if int(root.attrib.get('AnswerCount')) == 1:
                                    single_answer_questions += 1
                                continue # skip this row (line)

                            # add all the other attributes (without filtering)
                            for attr_id in Q_attr_to_keep:
                                row_dict[attr_id] = root.attrib.get(attr_id)
                                
                            questions_row_list.append(row_dict)

                        # keeping all answers, then we filter out those related
                        # to dropped questions (using 'filtered_ids.csv' file)
                        elif root.attrib.get('PostTypeId') == '2': # answer
                            row_dict = {}
                            for attr_id in A_attr_to_keep:
                                row_dict[attr_id] = root.attrib.get(attr_id)
                                
                            answers_row_list.append(row_dict)
                            
                        else:
                            pass # exclude entries neither questions nor answers
        
                lines_count += linesxchunk
                if lines_count / num_lines * 100 > prev_prog + 2:
                    if progress_out: 
                        print(f'...progress: {(lines_count / num_lines * 100):.2f}%...'.center(width), flush=True)
                    prev_prog = lines_count / num_lines * 100

        # convert the filtered data to Pandas dataframes
        df_questions = pd.DataFrame(questions_row_list)
        df_answers   = pd.DataFrame(answers_row_list)

        # keep only answers to the questions we have selected
        df_answers = df_answers[df_answers['ParentId'].isin(df_questions['Id'])]

        df_questions['timestamp'] = pd.to_datetime(df_questions['CreationDate'], format='%Y-%m-%dT%H:%M:%S.%f')
        df_answers['timestamp'] = pd.to_datetime(df_answers['CreationDate'], format='%Y-%m-%dT%H:%M:%S.%f')

        df_questions['Body'] = df_questions['Body'].apply(markdown_to_text) # convert HTML to plain text

        # saving the Pandas dataframes in 'csv' or 'pkl' format
        if file_format == 'csv':
            df_questions.to_csv(out_q, index=False, quoting=1, encoding='utf-8', errors='replace') # quote the fields
            df_answers.to_csv(out_a, index=False, encoding='utf-8', errors='replace')
        else:
            df_questions.to_pickle(out_q)
            df_answers.to_pickle(out_a)
        
        if th_num_answ > 0:
            print(f'!!! removed questions with < {th_num_answ} answers: {single_answer_questions}/{num_lines} !!!'.center(width), flush=True)
        
        print('... saved Q&A ...'.center(width), '\n', flush=True)

        ids_questions = df_questions['Id'].unique().tolist()
        ids_answers = df_answers['Id'].unique().tolist()
        all_ids = ids_questions + ids_answers # to then filter the votes
        df_all_ids = pd.DataFrame(all_ids, columns=['Id'])
        df_all_ids.to_csv(out_ids, index=False, encoding='utf-8', errors='replace')

        del df_questions # free memory
        del df_answers

    else: 
        print('... Q&A data already loaded ...'.center(width), '\n', flush=True)


    # ----------------------------- votes --------------------------------------

    out_v = os.path.join(out_folder, 'df_votes.' + file_format)

    if not os.path.exists(out_v):

        # use the csv of the ids to keep only the relevant records
        ids = set(pd.read_csv(out_ids)['Id'].tolist())

        votes_row_list = [] # list of dictionaries
        attr_to_keep = ['PostId', 'VoteTypeId', 'CreationDate']
        
        print('... saving votes to file ...'.center(width), flush=True)

        with open(out_votes, 'r', encoding='utf-8') as file:
            for _ in range(2): # (dataset-specific) skip first two lines
                next(file)

            for chunk in read_in_chunks(file, linesxchunk): # read chunks 
                for line in chunk:
                    if line.strip() != '</votes>': # signals end of file
                        root = ET.fromstring(line) # parse XML line
                        
                        if int(root.attrib.get('PostId')) in ids: # NOTE: inefficient
                            row_dict = {}
                            for attr_id in attr_to_keep:
                                row_dict[attr_id] = root.attrib.get(attr_id)
                            
                            votes_row_list.append(row_dict)

        df_votes = pd.DataFrame(votes_row_list)
        df_votes['timestamp'] = pd.to_datetime(df_votes['CreationDate'],
                                               format='%Y-%m-%dT%H:%M:%S.%f')
        
        if file_format == 'csv':    
            df_votes.to_csv(out_v, index=False, encoding='utf-8', errors='replace')
        else:
            df_votes.to_pickle(out_v)

        print('... saved votes ...'.center(width), '\n', flush=True)

        del df_votes # free memory

    else: 
        print('... votes data already loaded ...'.center(width), '\n', flush=True)

    
    # ----------------------------- badges -------------------------------------

    out_b = os.path.join(out_folder, 'df_badges.' + file_format)

    if not os.path.exists(out_b):

        badges_row_list = [] # list of dictionaries
        attr_to_keep = ['UserId', 'Name', 'Class']

        print('... saving badges to file ...'.center(width), flush=True)

        num_lines = sum(1 for _ in open(out_badge, 'rb'))
        lines_count, prev_prog = 0, -1 # textual progression (for server)

        with open(out_badge, 'r', encoding='utf-8') as file:
            for _ in range(2): # (dataset-specific) skip first two lines
                next(file)

            for chunk in read_in_chunks(file, linesxchunk): # read chunks 
                for line in chunk:
                    if line.strip() != '</badges>': # signals end of file
                        root = ET.fromstring(line)  # parse XML line
                        row_dict = {}
                        for attr_id in attr_to_keep:
                            row_dict[attr_id] = root.attrib.get(attr_id)
                            
                        badges_row_list.append(row_dict)

                lines_count += linesxchunk
                if lines_count / num_lines * 100 > prev_prog + 2:
                    if progress_out: 
                        print(f'...progress: {(lines_count / num_lines * 100):.2f}%...'.center(width), flush=True)
                    prev_prog = lines_count / num_lines * 100

        df_badges = pd.DataFrame(badges_row_list)

        if file_format == 'csv':
            df_badges.to_csv(out_b, index=False, encoding='utf-8', errors='replace')
        else:
            df_badges.to_pickle(out_b)
        
        print('... saved badges ...'.center(width), '\n', flush=True)

    else: 
        print('... badges data already loaded ...'.center(width), '\n', flush=True)
        
    # TODO: parse the tags xml, to have info on the questions discussed

    print('*** DATA LOADING COMPLETED ***'.center(width, '*'), '\n', flush=True)
    