# Perform the desired filtering on the dataset. In addition, we add the 'dif-
# ficulty' (Easy, Medium, Hard) of the question using the library Ollama to 
# perform zero-shot learning using pre-trained model (i.e. llama) on the basis
# of the 'Body' and 'Title' of the question.
# Note that the votes timestamp have a granularity of a day.

import os
import re
import shutil
import pandas as pd
from collections import defaultdict, Counter
import ollama

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


VERBOSE = True
MAX_WORDS_FOWARD = 20
ADD_DIFFICULTY = False


# Returns the difficulty level of a question (provided by the user as a title and 
# body) as a string ("Easy", "Medium", or "Hard").
# It is possible to select another model (if previously downloaded) and an
# alternative prompt, which consider a different point of view for the difficulty.
def zero_shot_classification(title_p, body_p, model_p='gemma2:2b', template_p='a'):

    # default prompt (difficult from a STEM student perspective)
    template = f"""
    Can you classify the following question as "Easy", "Medium", or "Hard", assuming
    a STEM student is asking the question? The title of the question is "{title_p}" 
    and the body of the question is "{body_p}". Please, just give a 1 word answer,
    corresponding to the difficulty level of the question. Do not provide an explanation.
    """
    candiates = ['easy', 'medium', 'hard']
    pattern = r'[^A-Za-z]'

    if template_p == 'b':
        template = f"""
        Can you classify the following question as "Easy", "Medium", or "Hard", assuming
        an individual with a background in the question's field is asking the question?
        The title of the question is "{title_p}" and the body of the question is "{body_p}".
        Please, just give a 1 word answer, corresponding to the difficulty level of the 
        question. Do not provide an explanation.
        """
    
    if template_p == 'c':
        template = f"""
        Can you classify the following question with a score from 1 to 5 according to 
        its difficulty (5 corresponds to most difficult), assuming
        an individual with a background in the question's field is asking the question?
        The title of the question is "{title_p}" and the body of the question is "{body_p}".
        Please, just give a number as the answer, corresponding to the difficulty level  
        of the question. Do not provide an explanation.
        """
        candiates = ['1', '2', '3', '4', '5']
        pattern = r'[^0-9]'

    output   = ollama.chat(model=model_p, messages=[{'role': 'user', 'content': template,},])
    response = output['message']['content']

    # remove newline characters and remove leading and trailing whitespaces
    cleaned_response = response.replace('\n', ' ').strip()

    # chck that the response is actually only one word, if not try to take the first word
    if len(cleaned_response.split()) > 1:
        if VERBOSE:
            print(f'Longer response: {cleaned_response}', '\n')

        # keep only letters (e.g. remove '**')
        candidate_response = re.sub(pattern, '', cleaned_response.split()[0])
        
        if candidate_response.lower() in candiates:
            return candidate_response

        # check if the classification is in the first sentence of the response
        else:
            for word in cleaned_response.split()[:MAX_WORDS_FOWARD]:
                candidate_response = re.sub(pattern, '', word)

                if candidate_response.lower() in candiates:
                    return candidate_response
                else:
                    if VERBOSE:
                        print(f'FATAL ERROR: No valid response found in:\n{response}')
                    with open('errors.txt', 'a', encoding='utf-8') as file: # append
                        file.write(response + '\t--------------\n','errors.txt')
                    return 'None'

    else:
        return re.sub(pattern, '', cleaned_response.split()[0])
    

def save_string_to_file(string, filename):
    with open(filename, 'a', encoding='utf-8') as file: # append
        file.write(string)


# saves a dictionary (to be used to save the info for the histograms)
# beware may create very large files if the number of unique values is high
def save_counter_data(counter, column_name, file_name):
    with open(file_name, 'w') as f:
        f.write(f"Counter data for {column_name}:\n")
        for key, value in counter.items():
            f.write(f"{key}: {value}\n")


# checks if a folder (or tree) exists, if not creates it, if it does it wipes it
def check_and_wipe(parent_name, dataset_p):
    if not os.path.isdir(os.path.join(parent_name, f'{dataset_p}')):
        os.makedirs(os.path.join(parent_name, f'{dataset_p}'))

    else: # wipe the previous info
        for filename in os.listdir(os.path.join(parent_name, f'{dataset_p}')):
            file_path = os.path.join(os.path.join(parent_name, f'{dataset_p}'), filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f'Failed to delete {file_path}. Reason: {e}')



# ------------------------------------ Main ------------------------------------

if __name__ == '__main__':

    dataset     = str(os.environ.get('DATASET', 'math'))
    file_format = str(os.environ.get('OUT_FORMAT', 'csv'))

    filter_by_viewcount = int(os.environ.get('FILTER_BY_VIEWCOUNT', 0))
    width = 80
    
    # read the tabular loaded data (either in csv of pkl format)
    try:
        base_path = os.path.join('.', 'data', file_format, dataset)

        if file_format == 'csv':
            df_questions = pd.read_csv(os.path.join(base_path, f'df_questions.{file_format}'))
            df_answers   = pd.read_csv(os.path.join(base_path, f'df_answers.{file_format}'))
            df_votes     = pd.read_csv(os.path.join(base_path, f'df_votes.{file_format}'))
            df_badges    = pd.read_csv(os.path.join(base_path, f'df_badges.{file_format}'))

        elif file_format == 'pkl':
            df_questions = pd.read_pickle(os.path.join(base_path, f'df_questions.{file_format}'))
            df_answers   = pd.read_pickle(os.path.join(base_path, f'df_answers.{file_format}'))
            df_votes     = pd.read_pickle(os.path.join(base_path, f'df_votes.{file_format}'))
            df_badges    = pd.read_pickle(os.path.join(base_path, f'df_badges.{file_format}'))

        else:
            raise 'ERROR: invalid file format'

    except FileNotFoundError: 
        raise 'ERROR: dataframe files missing, load the data first [load_dataset.py]'  
    

    # ------------------------- Filtering by viewcount -------------------------

    print('*** STARTED PREPROCESSING ***'.center(width, '*'), '\n', flush=True)
    print(f'Dataset: {dataset}, filtered by viewcount: {filter_by_viewcount}'.center(width, '-'),
           '\n', flush=True)
    
    pre = len(df_questions)
    df_questions = df_questions[df_questions['ViewCount'] >= filter_by_viewcount]

    print(f'Filtered [viewcount] Q df: {len(df_questions)}/{pre}'.center(width), '\n', flush=True)
    
    df_answers = df_answers[df_answers['ParentId'].isin(df_questions['Id'])]


    # ------------------ Collect overall statics of the data -------------------

    check_and_wipe('info', dataset)
    check_and_wipe('distrib', dataset)

    num_questions     = len(df_questions)
    nunique_questions = df_questions['Id'].nunique()
    if num_questions != nunique_questions:
        save_string_to_file(
            f'len() and nunique() differ for questions: {num_questions} vs {nunique_questions}.\n',
            os.path.join('info', f'{dataset}', 'warnings.txt')
        )

    num_answers     = len(df_answers)
    nunique_answers = df_answers['Id'].nunique()
    if num_answers != nunique_answers:
        save_string_to_file(
            f'len() and nunique() differ for answers: {num_answers} vs {nunique_answers}.\n',
            os.path.join('info', f'{dataset}', 'warnings.txt')
        )

    # determine also the number of users over the platform
    num_users = df_badges['UserId'].nunique()
    q_with_an_answ = len(df_questions[df_questions['AnswerCount']>0])

    save_string_to_file(
        f'In {dataset} there are {num_questions} questions and {num_answers} answers, with {num_users} unique users.\n' + 
        f'The questions with at least one answer are: {q_with_an_answ}/{num_questions}' + '\n',
        os.path.join('info', f'{dataset}', 'dataset.txt')
        )
    

    # save some info regarding the data distributions (Question/Answer score, ViewCount,
    # number of answers per question
    print('Collecting dataset metrics, may take a while'.center(width, '.'), '\n', flush=True)

    score_counter = Counter(df_questions['Score'].dropna())
    save_counter_data(score_counter, 'Score', 
                      os.path.join('distrib', f'{dataset}', 'questions_score.txt'))
    print('1/5 saved'.center(width), flush=True)

    viewcount_counter = Counter(df_questions['ViewCount'].dropna())
    save_counter_data(viewcount_counter, 'ViewCount', 
                      os.path.join('distrib', f'{dataset}', 'questions_views.txt'))
    print('2/5 saved'.center(width), flush=True)

    answer_count_counter = Counter(df_questions['AnswerCount'].dropna())
    save_counter_data(answer_count_counter, 'AnswerCount', 
                      os.path.join('distrib', f'{dataset}', 'questions_answers.txt'))
    print('3/5 saved'.center(width), flush=True)

    score_counter = Counter(df_answers['Score'].dropna())
    save_counter_data(score_counter, 'Score', 
                      os.path.join('distrib', f'{dataset}', 'answers_score.txt'))
    print('4/5 saved'.center(width), flush=True)

    viewcount_counter = Counter(df_answers['ViewCount'].dropna())
    save_counter_data(viewcount_counter, 'ViewCount', 
                      os.path.join('distrib', f'{dataset}', 'answers_views.txt'))
    print('5/5 saved'.center(width), '\n', flush=True)

    # TODO: save the tag distribution [need to load 'Tags.xml' file]

    # NOTE: additional metrics that may be collected: 1. number of up/down votes,
    #       2. number of answers/questions for each user, 3. reputation score


    # -------------------------- Best Badge for Users --------------------------
    
    badges_pre_filter = len(df_badges) # col: ['UserId', 'Name', 'Class']

    # best class for each user (Gold=1, Silver=2, Bronze=3)
    best_class = df_badges.groupby('UserId')[['Class']].min().reset_index()
    
    # keep only best class for each user (whatever the badge name) # TODO
    df_badges = df_badges.merge(
        best_class, 
        how='inner', 
        left_on='UserId', 
        right_on='UserId', 
        suffixes=(None, '_best')
    )

    df_badges = df_badges[df_badges['Class'] == df_badges['Class_best']].drop_duplicates('UserId').drop(columns=['Class_best'])
    
    diff = badges_pre_filter - len(df_badges)
    print(f'Dropped {diff}/{badges_pre_filter} badges from users with multiple ones'.center(width), 
          '\n', flush=True)

    
    # ---------------------- Zero-shot Classification -------------------------

    # add the difficulty level of the question
    model    = 'gemma2:2b'
    template = 'c'

    target_file = f'./data/{file_format}/{dataset}/df_qdiff_viewcount_{filter_by_viewcount}.{file_format}'

    if ADD_DIFFICULTY: 

        if os.path.exists(target_file):
            del df_questions
            df_questions = pd.read_csv(target_file) # load to avoid performing classification again
        
        else:
            print(f'Performing zero-shot classification with {model} - template {template}'.center(width), 
                '\n', flush=True)

            df_questions['difficulty'] = df_questions.apply(
                lambda x: zero_shot_classification(x['Title'], x['Body'],
                                                model_p=model, template_p=template),
                axis=1 # apply to each row
            )

            # save the dataframe with the difficulty score for each question
            df_questions.to_csv(target_file, index=False, encoding='utf-8', errors='ignore')


    # ------------------------- Get a Single DataFrame -------------------------

    # 'inner' merge keeps only rows with a matching value in both df (`left_on` `right_on`)

    # identifier of a post: 'Id' (question/answer), 'PostId' (votes)

    record_votes_pre_filter = len(df_votes)
    questions_pre_filter    = df_questions['Id'].nunique()
    answers_pre_filter      = df_answers['Id'].nunique()


    # merge 'vote' with 'answer' to have the 'ParentId' (in 'answer') in each vote 
    # record as this information is not provided in the XML vote's file. This df
    # contains the dynamics of the answers and the question they refer to.
    tmp = df_votes.merge(
        df_answers, 
        left_on='PostId', 
        right_on='Id', 
        how='inner', 
        suffixes=("_vote", "_answer") # suffixes for the conflicting columns
    )

    remaining_answers   = tmp['Id'].nunique()
    remaining_questions = tmp['ParentId'].nunique()

    print(f'After merging votes with answers:'.center(width, '.'), '\n', 
          f'{len(tmp)}/{record_votes_pre_filter} votes records'.center(width), '\n',
          f'{remaining_answers}/{answers_pre_filter} answers remaining'.center(width), '\n',
          f'{remaining_questions}/{questions_pre_filter} questions remaining'.center(width),
          '\n', flush=True)
    
    # NOTE: most of the 'missing answers' (the discrepancy in the info message above)
    #       is because (most) of those answers do NOT have a dynamics (up/down votes)
    # Zero-score answers are a good proxy for this (not the same thing because a zero
    # score may result from positive and negative votes that cancel each other out).
    df_zero_score_answers = df_answers[df_answers['Score'] == 0]

    # check how much of the zero-score answers are in the 'missing answers', for
    # those we can conlude that they do not have a dynamics
    df_missing_answers = df_answers[~df_answers['Id'].isin(tmp['Id'])] # what is not there

    set_missing    = set(df_missing_answers['Id'])
    set_zero_score = set(df_zero_score_answers['Id'])
    missing_no_dynamics = len(set_missing.intersection(set_zero_score))

    diff = answers_pre_filter-remaining_answers
    msg = f'For the missing answers, {missing_no_dynamics}/{diff} have score equal to 0.'
    print(msg.center(width), flush=True)
    save_string_to_file(msg + '\n', os.path.join('info', f'{dataset}', 'dataset.txt'))

    del df_zero_score_answers # free space
    del df_missing_answers


    # at the same time we also lose questions, this because some questions do not
    # have dynamics themselves, so potentially no record in votes
    missing_questions = set(df_questions['Id']).difference(set(tmp['ParentId']))

    df_questions_no_answer = df_questions[df_questions['AnswerCount'] == 0]
    missing_no_answers = missing_questions.intersection(set(df_questions_no_answer['Id']))
    
    msg = f'Of the missing questions {len(missing_no_answers)}/{len(missing_questions)} have no answer (no dynamics).'
    print(msg.center(width), flush=True)
    save_string_to_file(msg + '\n', os.path.join('info', f'{dataset}', 'dataset.txt'))

    # save some statistics of the no-answer questions (score and viewcount)
    no_answer_score_counter = Counter(df_questions_no_answer['Score'].dropna())
    save_counter_data(no_answer_score_counter, 'Score', 
                      os.path.join('distrib', f'{dataset}', 'questions_no_answers_score.txt'))
    
    no_answer_viewcount_counter = Counter(df_questions_no_answer['ViewCount'].dropna())
    save_counter_data(no_answer_viewcount_counter, 'ViewCount', 
                      os.path.join('distrib', f'{dataset}', 'questions_no_answers_views.txt'))
    

    set_questions_lost = missing_questions.difference(missing_no_answers)

    # why do we have additional excluded questions? let us look at the associated
    # answers, it happens we do not have the dynamics for those answers (in df_votes)
    # the merge drops questions for which we do not have matches in df_votes
    remaining = len(missing_questions) - len(missing_no_answers)

    answ_to_lost_questions = set(df_answers[df_answers['ParentId'].isin(set_questions_lost)]['Id'])
    df_answ_lost_questions = df_votes[df_votes['PostId'].isin(answ_to_lost_questions)]
    diff = remaining - df_answ_lost_questions['PostId'].nunique()
    msg = f'Of those remaining {diff}/{remaining}, we have no info of their answers in df_votes.'
    print(msg.center(width), flush=True)
    save_string_to_file(msg + '\n', os.path.join('info', f'{dataset}', 'dataset.txt'))

    del missing_questions
    del set_questions_lost
    del missing_no_answers
    del answ_to_lost_questions
    del df_questions_no_answer
    del df_answ_lost_questions


    # answers for which we do not have the 'votes' info but have score != 0, we do
    # not have their history, we shall filter them (also the corresponding questions)
    # filtering by their 'ParentId' (question) we can achieve this
    dyn_missing_answer = set_missing.difference(set_zero_score)
    save_string_to_file(str(dyn_missing_answer) + '\n', os.path.join('info', f'{dataset}', 'missing_answers.txt'))
    
    q_to_remove = set(df_answers[df_answers['Id'].isin(dyn_missing_answer)]['ParentId'])
    save_string_to_file(str(q_to_remove) + '\n', os.path.join('info', f'{dataset}', 'missing_questions.txt'))

    tmp = tmp[~tmp['ParentId'].isin(q_to_remove)] # remove both questions and answers

    remaining_answers   = tmp['Id'].nunique()
    remaining_questions = tmp['ParentId'].nunique()

    print('\n', f'After filtering no-votes questions/answers:'.center(width), '\n', 
          f'{len(tmp)}/{record_votes_pre_filter} votes records'.center(width), '\n',
          f'{remaining_answers}/{answers_pre_filter} answers remaining'.center(width),'\n',
          f'{remaining_questions}/{questions_pre_filter} questions remaining'.center(width),
          '\n', flush=True)

    del dyn_missing_answer
    del q_to_remove


    # recall that 'df_answer' and 'df_questions' come from the save XML file and
    # so have the same structure and headers. `tmp` contains all the 'df_answers'
    # info. We are merging a subset of 'df_questions' ('_question') with 'tmp' 
    # ('_answer'); we use the prefixes accordingly. Note also that thanks to the
    # above merge we will not have conflicting columns between questions and the
    # info in `tmp` from 'df_votes' as have already have been resolved with the
    # merge with 'df_answers' that have the same structure!

    sub_col_to_keep = ['Id', 'OwnerUserId', 'timestamp', 'AnswerCount'] # TODO: add 'Score' (question)
    tmp = df_questions[sub_col_to_keep].merge(
        tmp,
        left_on='Id', 
        right_on='ParentId', 
        how='inner', 
        suffixes=("_question", "_answer") # 'Id' is in both dataframes
    )
    # 'timestamp' in the resulting df will be that of the question

    remaining_questions = tmp['Id_question'].nunique()

    print(f'After merging with questions:'.center(width, '.'), '\n', 
          f'{len(tmp)}/{record_votes_pre_filter} votes records'.center(width), '\n',
          f'{remaining_questions}/{questions_pre_filter} questions remaining'.center(width),
          '\n', flush=True)


    # merge with badges to have 'Class' and 'Name' for each record, do a 'left' merge
    # so that to retain all the info from the 'tmp' dataframe
    df_all = tmp.merge(
        df_badges,
        left_on='OwnerUserId_question', 
        right_on='UserId', 
        how='left'
    ).drop(columns=['UserId']) # drop 'UserId', we can use 'OwnerUserId_question'

    remaining_questions = df_all['Id_question'].nunique()

    print(f'After merging with badges:'.center(width, '.'), '\n', 
          f'{len(df_all)}/{record_votes_pre_filter} votes records'.center(width), '\n',
          f'{remaining_questions}/{questions_pre_filter} questions remaining'.center(width),
          '\n', flush=True)


    # ----------- Filtering (dropna) and Preprocess Merged Dataset  ------------

    pre_filter = len(df_all)

    col_sub = [
            'PostId', 'VoteTypeId', 'timestamp', 'timestamp_vote', 'Id_answer',
            'PostTypeId', 'timestamp_answer', 'Score', 'ParentId'
        ] # 'OwnerUserId_question' may contain nan (but we do not use it for now)
    df_all = df_all.dropna(subset=col_sub)

    # do not cast all some may contain nan (e.g. 'Title')
    df_all = df_all.astype(
        {
            'PostId': 'int32',
            'VoteTypeId': 'int32',
            'Id_answer': 'int32',
            'PostTypeId': 'int32',
            'Score': 'int32',                    # final score of the question
            'timestamp': 'datetime64[ns]',       # question timestamp
            'timestamp_vote': 'datetime64[ns]',
            'timestamp_answer': 'datetime64[ns]',
            'ParentId': 'int32',
        }
    )

    diff = pre_filter-len(df_all)
    print(f'Dropped {diff}/{pre_filter} rows with None values'.center(width), flush=True)


    # -------------- Compute the Score Evolution for each Answer ---------------

    # sort once for all, (some) other operations maintain the order (e.g. groupby)
    # recall that we started from the answers so df_all contains basically all the
    # evolution of the answers (with additional columns fo questions info)
    df_all = df_all.sort_values('timestamp_vote') 

    vote_increment = defaultdict(int, {2: 1, 3: -1}) # 0 if the 'key' is not present
    
    df_all['increment'] = df_all['VoteTypeId'].map(vote_increment)                         # sequence of increments
    df_all['score_history'] = df_all.groupby(['ParentId', 'PostId'])['increment'].cumsum() # score evolution timeline


    # ------------ Colelct some additional metrics and save results ------------

    # determine also static metrics (e.g. time to first answer). Do after filtering:
    # fairer comparison with "time to emerge" or other metrics. Disregard 'CreationDate'
    # which are the non-converted timestamps, 'timestamp' is the one of the question, 
    # 'timestamp_vote' and 'timestamp_answer' are those of votes and answers respectively

    # time to first answer
    quickest_answer = df_all.groupby('ParentId')['timestamp_answer'].transform('min')
    df_all['time_to_first_answer'] = quickest_answer - df_all['timestamp']

    # avoid oversampling, keep one delta timestamp for each question (df_all has the entire dynamics)
    group = df_all.groupby('ParentId').first().reset_index()

    filename = os.path.join('distrib', f'{dataset}', 'time_to_first_answer.csv')
    group['time_to_first_answer'].to_csv(filename, index=False) # timedelta, saved as a str

    del group


    # time to best answer (time of posting of the best, not the time it became the best)
    best_score = df_all.groupby('ParentId')['Score'].transform('max')
    df_best_answer = df_all[df_all['Score'] == best_score]           # keep only best evolution
    group = df_best_answer.groupby('ParentId').first().reset_index() # avoid oversampling

    group['time_to_best_answer'] = group['timestamp_answer'] - group['timestamp']

    filename = os.path.join('distrib', f'{dataset}', 'time_to_best_answer.csv')
    group['time_to_best_answer'].to_csv(filename, index=False) # timedelta, saved as a str

    del group

    print('*** ENDED PREPROCESSING ***'.center(width, '*'), '\n', flush=True)


    # save the overall dataframe (csv preferred format even if pkl preserves the types)
    if file_format == 'csv':
        df_all.to_csv(f'./data/{file_format}/{dataset}/df_all_viewcount_{filter_by_viewcount}.{file_format}', index=False, encoding='utf-8', errors='ignore')
    
    elif file_format == 'pkl':
        df_all.to_pickle(f'./data/{file_format}/{dataset}/df_all_viewcount_{filter_by_viewcount}.{file_format}')
    
    else:
        raise 'ERROR: invalid file format'
