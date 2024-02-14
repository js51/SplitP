import pytest
import os

@pytest.fixture(scope="session")
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope="session")
def test_cases_trees():
    cases = [
        {
            "newick_string": "(0:0.1,1:0.4);",
            "number_of_taxa": 2,
            "true_splits": [],
            "model": None,
        },
        {
            "newick_string": "((0:0.1,1:0.3),2:0.03);",
            "number_of_taxa": 3,
            "true_splits": [],
            "model": None,
        },
        {
            "newick_string": "((0,2),(1,3));",
            "number_of_taxa": 4,
            "true_splits": ["02|13"],
            "model": ("JC", 0.05),
        },
        {
            "number_of_taxa": 4,
            "true_splits": ["01|23"],
            "model": ("JC", 0.05),
        },
        {
            "newick_string": "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
            "number_of_taxa": 4,
            "JSON_tree": {
                "id": "F",
                "children": [
                    {"id": "A", "branch_length": 0.1},
                    {"id": "B", "branch_length": 0.2},
                    {
                        "id": "E",
                        "branch_length": 0.5,
                        "children": [
                            {"id": "C", "branch_length": 0.3},
                            {"id": "D", "branch_length": 0.4},
                        ],
                    },
                ],
            },
        },
        {
            "newick_string": "((((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6):0.4,(E:0.7,F:0.8):0.9):0.3,(G:0.9,H:0.9):0.6);",
            "number_of_taxa": 8,
            "true_splits": [
                'ABCDEF|GH', 
                'ABEFGH|CD', 
                'AB|CDEFGH', 
                'ABCDGH|EF',
                'ABCD|EFGH',
            ]
        },
        {
            "newick_string": 
                "(((Drosophila_melanogaster:1,Caenorhabditis_elegans:1):1,((Ciona_intestinalis:1,Ciona_savignyi:1):1,(((((((((Taeniopygia_guttata:1,Ficedula_albicollis:1):1,((Meleagris_gallopavo:1,Gallus_gallus:1):1,Anas_platyrhynchos:1):1):1,Pelodiscus_sinensis:1):1,Anolis_carolinensis:1):1,((((((Procavia_capensis:1,Loxodonta_africana:1):1,Echinops_telfairi:1):1,(Choloepus_hoffmanni:1,Dasypus_novemcinctus:1):1):1,((((Oryctolagus_cuniculus:1,Ochotona_princeps:1):1,((((Mus_musculus:1,Rattus_norvegicus:1):1,Dipodomys_ordii:1):1,Ictidomys_tridecemlineatus:1):1,Cavia_porcellus:1):1):1,(((Microcebus_murinus:1,Otolemur_garnettii:1):1,(((((Papio_anubis:1,Macaca_mulatta:1):1,Chlorocebus_sabaeus:1):1,((((Pan_troglodytes:1,Homo_sapiens:1):1,Gorilla_gorilla:1):1,Pongo_abelii:1):1,Nomascus_leucogenys:1):1):1,Callithrix_jacchus:1):1,Tarsius_syrichta:1):1):1,Tupaia_belangeri:1):1):1,((Sorex_araneus:1,Erinaceus_europaeus:1):1,(((Pteropus_vampyrus:1,Myotis_lucifugus:1):1,((((Mustela_putorius_furo:1,Ailuropoda_melanoleuca:1):1,Canis_familiaris:1):1,Felis_catus:1):1,Equus_caballus:1):1):1,((((Bos_taurus:1,Ovis_aries:1):1,Tursiops_truncatus:1):1,Vicugna_pacos:1):1,Sus_scrofa:1):1):1):1):1):1,((Macropus_eugenii:1,Sarcophilus_harrisii:1):1,Monodelphis_domestica:1):1):1,Ornithorhynchus_anatinus:1):1):1,Xenopus_tropicalis:1):1,Latimeria_chalumnae:1):1,(((Danio_rerio:1,Astyanax_mexicanus:1):1,(((Tetraodon_nigroviridis:1,Takifugu_rubripes:1):1,((((Poecilia_formosa:1,Xiphophorus_maculatus:1):1,Oryzias_latipes:1):1,Gasterosteus_aculeatus:1):1,Oreochromis_niloticus:1):1):1,Gadus_morhua:1):1):1,Lepisosteus_oculatus:1):1):1,Petromyzon_marinus:1):1):1):1,Saccharomyces_cerevisiae:1);",
            "number_of_taxa": 60,
            "JSON_tree": {
                "children": [
                    {
                        "id": "",
                        "children": [
                            {
                                "id": "",
                                "children": [
                                    {
                                        "id": "Drosophila_melanogaster",
                                        "branch_length": 1,
                                    },
                                    {
                                        "id": "Caenorhabditis_elegans",
                                        "branch_length": 1,
                                    },
                                ],
                                "branch_length": 1,
                            },
                            {
                                "id": "",
                                "children": [
                                    {
                                        "id": "",
                                        "children": [
                                            {
                                                "id": "Ciona_intestinalis",
                                                "branch_length": 1,
                                            },
                                            {
                                                "id": "Ciona_savignyi",
                                                "branch_length": 1,
                                            },
                                        ],
                                        "branch_length": 1,
                                    },
                                    {
                                        "id": "",
                                        "children": [
                                            {
                                                "id": "",
                                                "children": [
                                                    {
                                                        "id": "",
                                                        "children": [
                                                            {
                                                                "id": "",
                                                                "children": [
                                                                    {
                                                                        "id": "",
                                                                        "children": [
                                                                            {
                                                                                "id": "",
                                                                                "children": [
                                                                                    {
                                                                                        "id": "",
                                                                                        "children": [
                                                                                            {
                                                                                                "id": "",
                                                                                                "children": [
                                                                                                    {
                                                                                                        "id": "",
                                                                                                        "children": [
                                                                                                            {
                                                                                                                "id": "Taeniopygia_guttata",
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                            {
                                                                                                                "id": "Ficedula_albicollis",
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                        ],
                                                                                                        "branch_length": 1,
                                                                                                    },
                                                                                                    {
                                                                                                        "id": "",
                                                                                                        "children": [
                                                                                                            {
                                                                                                                "id": "",
                                                                                                                "children": [
                                                                                                                    {
                                                                                                                        "id": "Meleagris_gallopavo",
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                    {
                                                                                                                        "id": "Gallus_gallus",
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                ],
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                            {
                                                                                                                "id": "Anas_platyrhynchos",
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                        ],
                                                                                                        "branch_length": 1,
                                                                                                    },
                                                                                                ],
                                                                                                "branch_length": 1,
                                                                                            },
                                                                                            {
                                                                                                "id": "Pelodiscus_sinensis",
                                                                                                "branch_length": 1,
                                                                                            },
                                                                                        ],
                                                                                        "branch_length": 1,
                                                                                    },
                                                                                    {
                                                                                        "id": "Anolis_carolinensis",
                                                                                        "branch_length": 1,
                                                                                    },
                                                                                ],
                                                                                "branch_length": 1,
                                                                            },
                                                                            {
                                                                                "id": "",
                                                                                "children": [
                                                                                    {
                                                                                        "id": "",
                                                                                        "children": [
                                                                                            {
                                                                                                "id": "",
                                                                                                "children": [
                                                                                                    {
                                                                                                        "id": "",
                                                                                                        "children": [
                                                                                                            {
                                                                                                                "id": "",
                                                                                                                "children": [
                                                                                                                    {
                                                                                                                        "id": "",
                                                                                                                        "children": [
                                                                                                                            {
                                                                                                                                "id": "Procavia_capensis",
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                            {
                                                                                                                                "id": "Loxodonta_africana",
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                        ],
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                    {
                                                                                                                        "id": "Echinops_telfairi",
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                ],
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                            {
                                                                                                                "id": "",
                                                                                                                "children": [
                                                                                                                    {
                                                                                                                        "id": "Choloepus_hoffmanni",
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                    {
                                                                                                                        "id": "Dasypus_novemcinctus",
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                ],
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                        ],
                                                                                                        "branch_length": 1,
                                                                                                    },
                                                                                                    {
                                                                                                        "id": "",
                                                                                                        "children": [
                                                                                                            {
                                                                                                                "id": "",
                                                                                                                "children": [
                                                                                                                    {
                                                                                                                        "id": "",
                                                                                                                        "children": [
                                                                                                                            {
                                                                                                                                "id": "",
                                                                                                                                "children": [
                                                                                                                                    {
                                                                                                                                        "id": "Oryctolagus_cuniculus",
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                    {
                                                                                                                                        "id": "Ochotona_princeps",
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                ],
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                            {
                                                                                                                                "id": "",
                                                                                                                                "children": [
                                                                                                                                    {
                                                                                                                                        "id": "",
                                                                                                                                        "children": [
                                                                                                                                            {
                                                                                                                                                "id": "",
                                                                                                                                                "children": [
                                                                                                                                                    {
                                                                                                                                                        "id": "",
                                                                                                                                                        "children": [
                                                                                                                                                            {
                                                                                                                                                                "id": "Mus_musculus",
                                                                                                                                                                "branch_length": 1,
                                                                                                                                                            },
                                                                                                                                                            {
                                                                                                                                                                "id": "Rattus_norvegicus",
                                                                                                                                                                "branch_length": 1,
                                                                                                                                                            },
                                                                                                                                                        ],
                                                                                                                                                        "branch_length": 1,
                                                                                                                                                    },
                                                                                                                                                    {
                                                                                                                                                        "id": "Dipodomys_ordii",
                                                                                                                                                        "branch_length": 1,
                                                                                                                                                    },
                                                                                                                                                ],
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                            {
                                                                                                                                                "id": "Ictidomys_tridecemlineatus",
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                        ],
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                    {
                                                                                                                                        "id": "Cavia_porcellus",
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                ],
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                        ],
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                    {
                                                                                                                        "id": "",
                                                                                                                        "children": [
                                                                                                                            {
                                                                                                                                "id": "",
                                                                                                                                "children": [
                                                                                                                                    {
                                                                                                                                        "id": "",
                                                                                                                                        "children": [
                                                                                                                                            {
                                                                                                                                                "id": "Microcebus_murinus",
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                            {
                                                                                                                                                "id": "Otolemur_garnettii",
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                        ],
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                    {
                                                                                                                                        "id": "",
                                                                                                                                        "children": [
                                                                                                                                            {
                                                                                                                                                "id": "",
                                                                                                                                                "children": [
                                                                                                                                                    {
                                                                                                                                                        "id": "",
                                                                                                                                                        "children": [
                                                                                                                                                            {
                                                                                                                                                                "id": "",
                                                                                                                                                                "children": [
                                                                                                                                                                    {
                                                                                                                                                                        "id": "",
                                                                                                                                                                        "children": [
                                                                                                                                                                            {
                                                                                                                                                                                "id": "Papio_anubis",
                                                                                                                                                                                "branch_length": 1,
                                                                                                                                                                            },
                                                                                                                                                                            {
                                                                                                                                                                                "id": "Macaca_mulatta",
                                                                                                                                                                                "branch_length": 1,
                                                                                                                                                                            },
                                                                                                                                                                        ],
                                                                                                                                                                        "branch_length": 1,
                                                                                                                                                                    },
                                                                                                                                                                    {
                                                                                                                                                                        "id": "Chlorocebus_sabaeus",
                                                                                                                                                                        "branch_length": 1,
                                                                                                                                                                    },
                                                                                                                                                                ],
                                                                                                                                                                "branch_length": 1,
                                                                                                                                                            },
                                                                                                                                                            {
                                                                                                                                                                "id": "",
                                                                                                                                                                "children": [
                                                                                                                                                                    {
                                                                                                                                                                        "id": "",
                                                                                                                                                                        "children": [
                                                                                                                                                                            {
                                                                                                                                                                                "id": "",
                                                                                                                                                                                "children": [
                                                                                                                                                                                    {
                                                                                                                                                                                        "id": "",
                                                                                                                                                                                        "children": [
                                                                                                                                                                                            {
                                                                                                                                                                                                "id": "Pan_troglodytes",
                                                                                                                                                                                                "branch_length": 1,
                                                                                                                                                                                            },
                                                                                                                                                                                            {
                                                                                                                                                                                                "id": "Homo_sapiens",
                                                                                                                                                                                                "branch_length": 1,
                                                                                                                                                                                            },
                                                                                                                                                                                        ],
                                                                                                                                                                                        "branch_length": 1,
                                                                                                                                                                                    },
                                                                                                                                                                                    {
                                                                                                                                                                                        "id": "Gorilla_gorilla",
                                                                                                                                                                                        "branch_length": 1,
                                                                                                                                                                                    },
                                                                                                                                                                                ],
                                                                                                                                                                                "branch_length": 1,
                                                                                                                                                                            },
                                                                                                                                                                            {
                                                                                                                                                                                "id": "Pongo_abelii",
                                                                                                                                                                                "branch_length": 1,
                                                                                                                                                                            },
                                                                                                                                                                        ],
                                                                                                                                                                        "branch_length": 1,
                                                                                                                                                                    },
                                                                                                                                                                    {
                                                                                                                                                                        "id": "Nomascus_leucogenys",
                                                                                                                                                                        "branch_length": 1,
                                                                                                                                                                    },
                                                                                                                                                                ],
                                                                                                                                                                "branch_length": 1,
                                                                                                                                                            },
                                                                                                                                                        ],
                                                                                                                                                        "branch_length": 1,
                                                                                                                                                    },
                                                                                                                                                    {
                                                                                                                                                        "id": "Callithrix_jacchus",
                                                                                                                                                        "branch_length": 1,
                                                                                                                                                    },
                                                                                                                                                ],
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                            {
                                                                                                                                                "id": "Tarsius_syrichta",
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                        ],
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                ],
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                            {
                                                                                                                                "id": "Tupaia_belangeri",
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                        ],
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                ],
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                            {
                                                                                                                "id": "",
                                                                                                                "children": [
                                                                                                                    {
                                                                                                                        "id": "",
                                                                                                                        "children": [
                                                                                                                            {
                                                                                                                                "id": "Sorex_araneus",
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                            {
                                                                                                                                "id": "Erinaceus_europaeus",
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                        ],
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                    {
                                                                                                                        "id": "",
                                                                                                                        "children": [
                                                                                                                            {
                                                                                                                                "id": "",
                                                                                                                                "children": [
                                                                                                                                    {
                                                                                                                                        "id": "",
                                                                                                                                        "children": [
                                                                                                                                            {
                                                                                                                                                "id": "Pteropus_vampyrus",
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                            {
                                                                                                                                                "id": "Myotis_lucifugus",
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                        ],
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                    {
                                                                                                                                        "id": "",
                                                                                                                                        "children": [
                                                                                                                                            {
                                                                                                                                                "id": "",
                                                                                                                                                "children": [
                                                                                                                                                    {
                                                                                                                                                        "id": "",
                                                                                                                                                        "children": [
                                                                                                                                                            {
                                                                                                                                                                "id": "",
                                                                                                                                                                "children": [
                                                                                                                                                                    {
                                                                                                                                                                        "id": "Mustela_putorius_furo",
                                                                                                                                                                        "branch_length": 1,
                                                                                                                                                                    },
                                                                                                                                                                    {
                                                                                                                                                                        "id": "Ailuropoda_melanoleuca",
                                                                                                                                                                        "branch_length": 1,
                                                                                                                                                                    },
                                                                                                                                                                ],
                                                                                                                                                                "branch_length": 1,
                                                                                                                                                            },
                                                                                                                                                            {
                                                                                                                                                                "id": "Canis_familiaris",
                                                                                                                                                                "branch_length": 1,
                                                                                                                                                            },
                                                                                                                                                        ],
                                                                                                                                                        "branch_length": 1,
                                                                                                                                                    },
                                                                                                                                                    {
                                                                                                                                                        "id": "Felis_catus",
                                                                                                                                                        "branch_length": 1,
                                                                                                                                                    },
                                                                                                                                                ],
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                            {
                                                                                                                                                "id": "Equus_caballus",
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                        ],
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                ],
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                            {
                                                                                                                                "id": "",
                                                                                                                                "children": [
                                                                                                                                    {
                                                                                                                                        "id": "",
                                                                                                                                        "children": [
                                                                                                                                            {
                                                                                                                                                "id": "",
                                                                                                                                                "children": [
                                                                                                                                                    {
                                                                                                                                                        "id": "",
                                                                                                                                                        "children": [
                                                                                                                                                            {
                                                                                                                                                                "id": "Bos_taurus",
                                                                                                                                                                "branch_length": 1,
                                                                                                                                                            },
                                                                                                                                                            {
                                                                                                                                                                "id": "Ovis_aries",
                                                                                                                                                                "branch_length": 1,
                                                                                                                                                            },
                                                                                                                                                        ],
                                                                                                                                                        "branch_length": 1,
                                                                                                                                                    },
                                                                                                                                                    {
                                                                                                                                                        "id": "Tursiops_truncatus",
                                                                                                                                                        "branch_length": 1,
                                                                                                                                                    },
                                                                                                                                                ],
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                            {
                                                                                                                                                "id": "Vicugna_pacos",
                                                                                                                                                "branch_length": 1,
                                                                                                                                            },
                                                                                                                                        ],
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                    {
                                                                                                                                        "id": "Sus_scrofa",
                                                                                                                                        "branch_length": 1,
                                                                                                                                    },
                                                                                                                                ],
                                                                                                                                "branch_length": 1,
                                                                                                                            },
                                                                                                                        ],
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                ],
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                        ],
                                                                                                        "branch_length": 1,
                                                                                                    },
                                                                                                ],
                                                                                                "branch_length": 1,
                                                                                            },
                                                                                            {
                                                                                                "id": "",
                                                                                                "children": [
                                                                                                    {
                                                                                                        "id": "",
                                                                                                        "children": [
                                                                                                            {
                                                                                                                "id": "Macropus_eugenii",
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                            {
                                                                                                                "id": "Sarcophilus_harrisii",
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                        ],
                                                                                                        "branch_length": 1,
                                                                                                    },
                                                                                                    {
                                                                                                        "id": "Monodelphis_domestica",
                                                                                                        "branch_length": 1,
                                                                                                    },
                                                                                                ],
                                                                                                "branch_length": 1,
                                                                                            },
                                                                                        ],
                                                                                        "branch_length": 1,
                                                                                    },
                                                                                    {
                                                                                        "id": "Ornithorhynchus_anatinus",
                                                                                        "branch_length": 1,
                                                                                    },
                                                                                ],
                                                                                "branch_length": 1,
                                                                            },
                                                                        ],
                                                                        "branch_length": 1,
                                                                    },
                                                                    {
                                                                        "id": "Xenopus_tropicalis",
                                                                        "branch_length": 1,
                                                                    },
                                                                ],
                                                                "branch_length": 1,
                                                            },
                                                            {
                                                                "id": "Latimeria_chalumnae",
                                                                "branch_length": 1,
                                                            },
                                                        ],
                                                        "branch_length": 1,
                                                    },
                                                    {
                                                        "id": "",
                                                        "children": [
                                                            {
                                                                "id": "",
                                                                "children": [
                                                                    {
                                                                        "id": "",
                                                                        "children": [
                                                                            {
                                                                                "id": "Danio_rerio",
                                                                                "branch_length": 1,
                                                                            },
                                                                            {
                                                                                "id": "Astyanax_mexicanus",
                                                                                "branch_length": 1,
                                                                            },
                                                                        ],
                                                                        "branch_length": 1,
                                                                    },
                                                                    {
                                                                        "id": "",
                                                                        "children": [
                                                                            {
                                                                                "id": "",
                                                                                "children": [
                                                                                    {
                                                                                        "id": "",
                                                                                        "children": [
                                                                                            {
                                                                                                "id": "Tetraodon_nigroviridis",
                                                                                                "branch_length": 1,
                                                                                            },
                                                                                            {
                                                                                                "id": "Takifugu_rubripes",
                                                                                                "branch_length": 1,
                                                                                            },
                                                                                        ],
                                                                                        "branch_length": 1,
                                                                                    },
                                                                                    {
                                                                                        "id": "",
                                                                                        "children": [
                                                                                            {
                                                                                                "id": "",
                                                                                                "children": [
                                                                                                    {
                                                                                                        "id": "",
                                                                                                        "children": [
                                                                                                            {
                                                                                                                "id": "",
                                                                                                                "children": [
                                                                                                                    {
                                                                                                                        "id": "Poecilia_formosa",
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                    {
                                                                                                                        "id": "Xiphophorus_maculatus",
                                                                                                                        "branch_length": 1,
                                                                                                                    },
                                                                                                                ],
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                            {
                                                                                                                "id": "Oryzias_latipes",
                                                                                                                "branch_length": 1,
                                                                                                            },
                                                                                                        ],
                                                                                                        "branch_length": 1,
                                                                                                    },
                                                                                                    {
                                                                                                        "id": "Gasterosteus_aculeatus",
                                                                                                        "branch_length": 1,
                                                                                                    },
                                                                                                ],
                                                                                                "branch_length": 1,
                                                                                            },
                                                                                            {
                                                                                                "id": "Oreochromis_niloticus",
                                                                                                "branch_length": 1,
                                                                                            },
                                                                                        ],
                                                                                        "branch_length": 1,
                                                                                    },
                                                                                ],
                                                                                "branch_length": 1,
                                                                            },
                                                                            {
                                                                                "id": "Gadus_morhua",
                                                                                "branch_length": 1,
                                                                            },
                                                                        ],
                                                                        "branch_length": 1,
                                                                    },
                                                                ],
                                                                "branch_length": 1,
                                                            },
                                                            {
                                                                "id": "Lepisosteus_oculatus",
                                                                "branch_length": 1,
                                                            },
                                                        ],
                                                        "branch_length": 1,
                                                    },
                                                ],
                                                "branch_length": 1,
                                            },
                                            {
                                                "id": "Petromyzon_marinus",
                                                "branch_length": 1,
                                            },
                                        ],
                                        "branch_length": 1,
                                    },
                                ],
                                "branch_length": 1,
                            },
                        ],
                        "branch_length": 1,
                    },
                    {"id": "Saccharomyces_cerevisiae", "branch_length": 1},
                ],
                "id": "",
            },
        },
    ]
    return cases