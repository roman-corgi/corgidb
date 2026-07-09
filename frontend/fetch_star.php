<?php 
header("Content-Type:application/json");
include("config.php"); ?>
<?php

// Check inputs 
if (!isset($_GET['st_name'])){
    die("st_name not set");
}

$conn = new mysqli($servername, $username, $password, $dbname);
// Check connection
if ($conn->connect_error) {
    die("Connection failed: " . $conn->connect_error);
}

$st_name = $_GET['st_name'];

$alias_result = $conn->query("SELECT `main_id` from StarAliases where `alias` = '$st_name'");
$alias_row = $alias_result->fetch_row();
$main_id = $alias_row[0];

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
    from Stars where main_id = '$main_id'";

$result = $conn->query($sql);
$data = $result->fetch_all();

// Return JSON
echo json_encode($data);

$conn->close();
?>

