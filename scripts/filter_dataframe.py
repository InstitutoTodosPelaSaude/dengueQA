import pandas as pd
import time
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter dataframe based on user defined parameters provided in a \'config.tsv\' file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--config", required=True, help="Parameter file")
    args = parser.parse_args()

    config = args.config # Example: /Users/anderson/google_drive/ITpS/projetos_colaboracoes/phyloDF/data/genomics/config_filter.tsv

    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df

    params = load_table(config)
    params.fillna('', inplace=True)
    params = params.set_index('param')

    # load data
    input_file = params.loc['input', 'value']
    df = load_table(input_file)


    # filter rows
    def filter_df(df, criteria):
        print('\nFiltering rows...')
        # print(criteria)
        new_df = pd.DataFrame()
        include = {}
        for filter_value in criteria.split(','):
            filter_value = filter_value.strip()
            if not filter_value.startswith('~'):
                col, val = filter_value.split(':')[0], filter_value.split(':')[1]
                if val == '\'\'':
                    val = ''
                if col not in include:
                    include[col] = [val]
                else:
                    include[col].append(val)
        # print('Include:', include)
        for filter_col, filter_val in include.items():
            print('\t- Including only rows with \'' + filter_col + '\' = \'' + ', '.join(filter_val) + '\'')
            # print(new_df.size)
            if new_df.empty:
                df_filtered = df[df[filter_col].isin(filter_val)]
                new_df = new_df.append(df_filtered)
            else:
                new_df = new_df[new_df[filter_col].isin(filter_val)]
            # print(new_df)#.head())

        exclude = {}
        for filter_value in criteria.split(','):
            filter_value = filter_value.strip()
            if filter_value.startswith('~'):
                # print('\t- Excluding all rows with \'' + col + '\' = \'' + val + '\'')
                filter_value = filter_value[1:]
                col, val = filter_value.split(':')[0], filter_value.split(':')[1]
                if val == '\'\'':
                    val = ''
                if col not in exclude:
                    exclude[col] = [val]
                else:
                    exclude[col].append(val)
        # print('Exclude:', exclude)
        for filter_col, filter_val in exclude.items():
            print('\t- Excluding all rows with \'' + filter_col + '\' = \'' + ', '.join(filter_val) + '\'')
            if new_df.empty:
                df = df[~df[filter_col].isin(filter_val)]
                new_df = new_df.append(df)
            else:
                new_df = new_df[~new_df[filter_col].isin(filter_val)]
            # print(new_df)#.head())
        return new_df

    # load data
    filters = params.loc['filter', 'value']
    if filters not in ['', None]:
        df = filter_df(df, filters)

    # filter by date
    date_col = params.loc['date_column', 'value']
    if date_col not in ['', None]:
        print('\t- Removing rows with incomplete dates')
        # remove genomes with incomplete dates
        df = df[df[date_col].apply(lambda x: len(x.split('-')) == 3)]  # accept only full dates
        df = df[df[date_col].apply(lambda x: 'X' not in x)]  # exclude -XX-XX missing dates

        start_date = params.loc['start_date', 'value']
        end_date = params.loc['end_date', 'value']

        today = time.strftime('%Y-%m-%d', time.gmtime())
        df[date_col] = pd.to_datetime(df[date_col])  # converting to datetime format
        if start_date in ['', None]:
            start_date = df[date_col].min()
        if end_date in ['', None]:
            end_date = today

        print('\t- Adding rows with %s from %s to %s' % (date_col, start_date, end_date))

        # converting back to string
        df[date_col] = df[date_col].apply(lambda x: x.strftime('%Y-%m-%d'))

        # filter genomes based on sampling date
        def filter_bydate(df, date):
            df[date] = pd.to_datetime(df[date])  # converting to datetime format
            mask = (df[date] >= start_date) & (df[date] <= end_date)  # mask any lines with dates outside the start/end dates
            df = df.loc[mask]  # apply mask
            return df

        df = filter_bydate(df, date_col)


    # drop columns
    ignore_cols = params.loc['ignore_cols', 'value']
    if ignore_cols not in ['', None]:
        for col in ignore_cols.split(','):
            df = df.drop(columns=col.strip())

    # sort data
    sort_by = params.loc['sort_by', 'value']
    if sort_by not in ['', None]:
        df = df.sort_values(by=sort_by)


    # output combined dataframe
    output_file = params.loc['output', 'value']
    df.to_csv(output_file, sep='\t', index=False)
    print('\nData successfully filtered and saved in:\n%s\n' % output_file)
