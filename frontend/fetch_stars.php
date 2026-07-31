<?php
header("Content-Type:application/json");
include("config.php"); ?>
<?php

// Check inputs
if (!isset($_GET['st_names'])){
    die("st_names not set");
}

$conn = new mysqli($servername, $username, $password, $dbname);
// Check connection
if ($conn->connect_error) {
    die("Connection failed: " . $conn->connect_error);
}

$st_names = explode(',', $_GET['st_names']);

$main_ids = array();
foreach ($st_names as $st_name) {
    $st_name = trim($st_name);
    $alias_result = $conn->query("SELECT `main_id` from StarAliases where `alias` = '$st_name'");
    $alias_row = $alias_result->fetch_row();
    if ($alias_row) {
        $main_ids[] = $alias_row[0];
    }
}

if (empty($main_ids)) {
    echo json_encode(array());
    $conn->close();
    exit;
}

$quoted_ids = array_map(function ($id) {
    return "'" . $id . "'";
}, $main_ids);
$id_list = implode(',', $quoted_ids);

$sql = "SELECT
    `st_name`,
    `main_id`,
    `ra`,
    `dec`,
    `spectype`,
    `sy_vmag`,
    `sy_imag`,
    `sy_dist`,
    `sy_plx`,
    `sy_pmra`,
    `sy_pmdec`,
    `st_radv`
    from Stars where main_id in ($id_list)";

$result = $conn->query($sql);
$data = $result->fetch_all();

// Return JSON
echo json_encode($data);

$conn->close();
?>
