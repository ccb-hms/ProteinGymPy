import pandas 

def getBioPlex(cell_line, version):
    '''
    Load BioPlex interactions data.
    
    This function loads BioPlex PPI data for
    cell lines HEK293T and HCT116, note we
    only have version 1.0 for HCT116 cells.
    
    Parameters
    ----------
    cell_line : str
        Takes input ['293T','HCT116'].
    version : str
        Takes input ['3.0','1.0','2.0'].
    
    Returns
    -------
    Pandas DataFrame
        A dataframe with each row corresponding to a PPI interaction.
    
    Examples
    --------
    >>> bp_293t = getBioPlex('293T', '1.0')
    >>> bp_293t.head(1)
       GeneA   GeneB UniprotA UniprotB SymbolA SymbolB            pW       pNI      pInt
    0    100  728378   P00813   A5A3E0     ADA   POTEF  2.605947e-09  0.000333  0.999667
    '''
    if (f'{cell_line}.{version}' not in 
        ['293T.1.0','293T.2.0','293T.3.0','HCT116.1.0']):
        print('dataset not available for this Cell Line - Version')
        
    else:
        if f'{cell_line}.{version}' == '293T.1.0':
            file_ext = 'interactionList_v2'
        elif f'{cell_line}.{version}' == '293T.2.0':
            file_ext = 'interactionList_v4a'
        elif f'{cell_line}.{version}' == '293T.3.0':
            file_ext = '293T_Network_10K_Dec_2019'
        elif f'{cell_line}.{version}' == 'HCT116.1.0':
            file_ext = 'HCT116_Network_5.5K_Dec_2019'

        BioPlex_interactions_df = pd.read_csv(
                f"https://bioplex.hms.harvard.edu/data/BioPlex_{file_ext}.tsv", 
                sep = '\t')
        
    # if pulling 293T cell line version 1.0 or 2.0, change column names to 
    # standardize across datasets for input into other functions
    if (cell_line == '293T') and (version == '1.0'):
        BioPlex_interactions_df.rename({'Gene A':'GeneA','Gene B':'GeneB',
            'Uniprot A':'UniprotA','Uniprot B':'UniprotB',
            'Symbol A':'SymbolA','Symbol B':'SymbolB','p(Wrong)':'pW',
            'p(No Interaction)':'pNI','p(Interaction)':'pInt'}, 
            axis = 1, inplace = True)

    if (cell_line == '293T') and (version == '2.0'):
        BioPlex_interactions_df.rename({'p(Wrong)':'pW',
            'p(No Interaction)':'pNI','p(Interaction)':'pInt'}, 
            axis = 1, inplace = True)
    
    return BioPlex_interactions_df


def getBioPlex2(cell_line, version):
    '''
    Load BioPlex interactions data.
    
    This function loads BioPlex PPI data for
    cell lines HEK293T and HCT116, note we
    only have version 1.0 for HCT116 cells.
    
    Parameters
    ----------
    cell_line : str
        Takes input ['293T','HCT116'].
    version : str
        Takes input ['3.0','1.0','2.0'].
    
    Returns
    -------
    Pandas DataFrame
        A dataframe with each row corresponding to a PPI interaction.
    
    Examples
    --------
    >>> bp_293t = getBioPlex('293T', '1.0')
    >>> bp_293t.head(1)
       GeneA   GeneB UniprotA UniprotB SymbolA SymbolB            pW       pNI      pInt
    0    100  728378   P00813   A5A3E0     ADA   POTEF  2.605947e-09  0.000333  0.999667
    '''
    if (f'{cell_line}.{version}' not in 
        ['293T.1.0','293T.2.0','293T.3.0','HCT116.1.0']):
        print('dataset not available for this Cell Line - Version')
        
    else:
        if f'{cell_line}.{version}' == '293T.1.0':
            file_ext = 'interactionList_v2'
        elif f'{cell_line}.{version}' == '293T.2.0':
            file_ext = 'interactionList_v4a'
        elif f'{cell_line}.{version}' == '293T.3.0':
            file_ext = '293T_Network_10K_Dec_2019'
        elif f'{cell_line}.{version}' == 'HCT116.1.0':
            file_ext = 'HCT116_Network_5.5K_Dec_2019'

        BioPlex_interactions_df = pd.read_csv(
                f"https://bioplex.hms.harvard.edu/data/BioPlex_{file_ext}.tsv", 
                sep = '\t')
        
    # if pulling 293T cell line version 1.0 or 2.0, change column names to 
    # standardize across datasets for input into other functions
    if (cell_line == '293T') and (version == '1.0'):
        BioPlex_interactions_df.rename({'Gene A':'GeneA','Gene B':'GeneB',
            'Uniprot A':'UniprotA','Uniprot B':'UniprotB',
            'Symbol A':'SymbolA','Symbol B':'SymbolB','p(Wrong)':'pW',
            'p(No Interaction)':'pNI','p(Interaction)':'pInt'}, 
            axis = 1, inplace = True)

    if (cell_line == '293T') and (version == '2.0'):
        BioPlex_interactions_df.rename({'p(Wrong)':'pW',
            'p(No Interaction)':'pNI','p(Interaction)':'pInt'}, 
            axis = 1, inplace = True)
    
    return BioPlex_interactions_df

