"""
Pytest configuration and fixtures for SEgene Analyzer tests
"""
import pytest
import tempfile
import shutil
import os
import pandas as pd

@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests"""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def sample_bed_file(temp_dir):
    """Create a sample BED file for testing"""
    bed_file = os.path.join(temp_dir, "test.bed")
    with open(bed_file, 'w') as f:
        f.write("chr1\t1000\t2000\tregion1\t100\n")
        f.write("chr2\t3000\t4000\tregion2\t200\n")
        f.write("chr7\t27000000\t28000000\tMYC_region\t500\n")
        f.write("chr8\t128000000\t129000000\tMYC_enhancer\t300\n")
    return bed_file

@pytest.fixture
def sample_info_file(temp_dir):
    """Create a sample information file for testing"""
    sample_file = os.path.join(temp_dir, "sample_info.txt")
    with open(sample_file, 'w') as f:
        f.write("sample_id\ttissue\tcell_type\tbiosample_type\n")
        f.write("sample1\tBrain\tNeuron\tPrimary\n")
        f.write("sample2\tLiver\tHepatocyte\tCulture\n")
        f.write("sample3\tHeart\tCardiomyocyte\tPrimary\n")
        f.write("sample4\tBrain\tAstrocyte\tPrimary\n")
        f.write("sample5\tLung\tPneumocyte\tCulture\n")
    return sample_file

@pytest.fixture
def regions_tsv_file(temp_dir):
    """Create a regions TSV file for testing"""
    regions_file = os.path.join(temp_dir, "regions.tsv")
    with open(regions_file, 'w') as f:
        f.write("chromosome\tstart\tend\tname\tdescription\n")
        f.write("chr1\t1500\t2500\ttest_region1\tTest region 1\n")
        f.write("chr2\t3500\t4500\ttest_region2\tTest region 2\n")
        f.write("chr7\t27500000\t27600000\tMYC_test\tMYC test region\n")
    return regions_file

@pytest.fixture
def regions_bed_file(temp_dir):
    """Create a regions BED file for testing"""
    regions_file = os.path.join(temp_dir, "regions.bed")
    with open(regions_file, 'w') as f:
        f.write("chr1\t1500\t2500\ttest_region1\n")
        f.write("chr2\t3500\t4500\ttest_region2\n")
        f.write("chr7\t27500000\t27600000\tMYC_test\n")
    return regions_file

@pytest.fixture
def sample_bed_dataframe():
    """Create a sample BED DataFrame for testing"""
    data = {
        'chr': ['chr1', 'chr2', 'chr7', 'chr8'],
        'start': [1000, 3000, 27000000, 128000000],
        'end': [2000, 4000, 28000000, 129000000],
        'name': ['region1', 'region2', 'MYC_region', 'MYC_enhancer'],
        'score': [100, 200, 500, 300]
    }
    return pd.DataFrame(data)

@pytest.fixture
def sample_info_dataframe():
    """Create a sample information DataFrame for testing"""
    data = {
        'sample_id': ['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
        'tissue': ['Brain', 'Liver', 'Heart', 'Brain', 'Lung'],
        'cell_type': ['Neuron', 'Hepatocyte', 'Cardiomyocyte', 'Astrocyte', 'Pneumocyte'],
        'biosample_type': ['Primary', 'Culture', 'Primary', 'Primary', 'Culture']
    }
    return pd.DataFrame(data)

@pytest.fixture
def mock_tissue_distribution():
    """Create mock tissue distribution data for testing"""
    return {
        'Brain': 25,
        'Liver': 15,
        'Heart': 10,
        'Lung': 12,
        'Kidney': 8,
        'Muscle': 6,
        'Skin': 5,
        'Others': 19
    }

@pytest.fixture
def mock_cell_id_counts():
    """Create mock cell ID counts for testing"""
    return {
        'ENCSR000EUA': 5,
        'ENCSR000EUB': 3,
        'ENCSR000EUC': 7,
        'ENCSR000EUD': 2,
        'ENCSR000EUE': 4
    }

@pytest.fixture
def output_dir(temp_dir):
    """Create an output directory for tests"""
    output_dir = os.path.join(temp_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    return output_dir