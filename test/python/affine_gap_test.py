import duckdb

def test_affine_gap():
    conn = duckdb.connect('');
    conn.execute("SELECT affine_gap('Sam') as value;");
    res = conn.fetchall()
    assert(res[0][0] == "Affine_gap Sam 🐥");