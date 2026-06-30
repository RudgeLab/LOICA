"""Regression tests for LOICA genetic network SBOL3 export."""

import pytest
import sbol3

from loica.geneproduct import Regulator, Reporter
from loica.genetic_network import GeneticNetwork
from loica.operators.hill1 import Hill1
from loica.operators.hill2 import Hill2
from loica.operators.receiver import Receiver
from loica.operators.source import Source
from loica.supplement import Supplement


@pytest.fixture(autouse=True)
def sbol_namespace():
    sbol3.set_namespace("https://example.org/loica-test")


def dna_part(display_id):
    component = sbol3.Component(display_id, sbol3.SBO_DNA)
    component.roles.append(sbol3.SO_ENGINEERED_REGION)
    return component


def regulator(name):
    return Regulator(name, sbol_comp=dna_part(f"{name}_coding_sequence"))


def reporter(name):
    return Reporter(name, sbol_comp=dna_part(f"{name}_coding_sequence"))


def supplement(name):
    return Supplement(name, sbol_comp=sbol3.Component(f"{name}_small_molecule", sbol3.SBO_SIMPLE_CHEMICAL))


def assert_exports_to_sbol(tmp_path, network, expected_min_objects):
    doc = network.to_sbol()

    assert isinstance(doc, sbol3.Document)
    assert len(doc.objects) >= expected_min_objects
    assert doc.find("geneticnetwork") is not None

    identities = [obj.identity for obj in doc.objects]
    report_sbol3 = doc.validate()
    assert len(report_sbol3) == 0, f"SBOL3 validation issues for {identities}: {report_sbol3}"

    output_file = tmp_path / "network.nt"
    doc.write(str(output_file), sbol3.SORTED_NTRIPLES)

    assert output_file.exists()
    assert output_file.stat().st_size > 0
    assert "geneticnetwork" in output_file.read_text()


@pytest.mark.parametrize(
    ("network", "expected_min_objects"),
    [
        pytest.param(
            GeneticNetwork(),
            1,
            id="empty-network",
        ),
    ],
)
def test_empty_network_exports(tmp_path, network, expected_min_objects):
    assert_exports_to_sbol(tmp_path, network, expected_min_objects)


def test_source_network_exports_to_sbol(tmp_path):
    gfp = reporter("gfp_source")
    source = Source(gfp, rate=10.0, name="source_promoter", sbol_comp=dna_part("source_promoter"))
    network = GeneticNetwork()
    network.add_reporter(gfp)
    network.add_operator(source)

    assert_exports_to_sbol(tmp_path, network, expected_min_objects=4)


def test_buffer_hill1_network_exports_to_sbol(tmp_path):
    input_reg = regulator("buffer_input")
    output_reg = regulator("buffer_output")
    buffer = Hill1(input_reg, output_reg, alpha=[1.0, 20.0], K=5.0, n=2.0, name="buffer", sbol_comp=dna_part("buffer_promoter"))
    network = GeneticNetwork()
    network.add_regulator([input_reg, output_reg])
    network.add_operator(buffer)

    assert_exports_to_sbol(tmp_path, network, expected_min_objects=5)


def test_not_hill1_network_exports_to_sbol(tmp_path):
    input_reg = regulator("not_input")
    output_reg = regulator("not_output")
    not_gate = Hill1(input_reg, output_reg, alpha=[20.0, 1.0], K=5.0, n=2.0, name="not_gate", sbol_comp=dna_part("not_promoter"))
    network = GeneticNetwork()
    network.add_regulator([input_reg, output_reg])
    network.add_operator(not_gate)

    assert_exports_to_sbol(tmp_path, network, expected_min_objects=5)


def test_receiver_network_exports_to_sbol(tmp_path):
    arabinose = supplement("arabinose")
    output_reg = regulator("receiver_output")
    receiver = Receiver(arabinose, output_reg, alpha=[1.0, 25.0], K=3.0, n=2.0, name="receiver", sbol_comp=dna_part("receiver_promoter"))
    network = GeneticNetwork()
    network.add_regulator(output_reg)
    network.add_operator(receiver)

    assert_exports_to_sbol(tmp_path, network, expected_min_objects=5)


def test_nor_hill2_network_exports_to_sbol(tmp_path):
    input_a = regulator("nor_input_a")
    input_b = regulator("nor_input_b")
    output_reg = regulator("nor_output")
    nor_gate = Hill2([input_a, input_b], output_reg, alpha=[20.0, 1.0, 1.0, 0.2], K=[5.0, 7.0], n=[2.0, 2.0], name="nor_gate", sbol_comp=dna_part("nor_promoter"))
    network = GeneticNetwork()
    network.add_regulator([input_a, input_b, output_reg])
    network.add_operator(nor_gate)

    assert_exports_to_sbol(tmp_path, network, expected_min_objects=6)


def test_repressilator_network_exports_to_sbol(tmp_path):
    lac_i = regulator("lacI")
    tet_r = regulator("tetR")
    c_i = regulator("cI")
    network = GeneticNetwork()
    network.add_regulator([lac_i, tet_r, c_i])
    network.add_operator([
        Hill1(c_i, lac_i, alpha=[20.0, 1.0], K=5.0, n=2.0, name="pLac", sbol_comp=dna_part("pLac")),
        Hill1(lac_i, tet_r, alpha=[20.0, 1.0], K=5.0, n=2.0, name="pTet", sbol_comp=dna_part("pTet")),
        Hill1(tet_r, c_i, alpha=[20.0, 1.0], K=5.0, n=2.0, name="pLambda", sbol_comp=dna_part("pLambda")),
    ])

    assert_exports_to_sbol(tmp_path, network, expected_min_objects=10)
