<?php 
header("Content-Type:application/json");
include("config.php"); ?>
<?php

// Check inputs
if (isset($_GET['st_names'])) {
    $st_names = explode(',', $_GET['st_names']);
} elseif (isset($_GET['st_name'])) {
    $st_names = array($_GET['st_name']);
} else {
    die("st_name or st_names not set");
}

$conn = new mysqli($servername, $username, $password, $dbname);
// Check connection
if ($conn->connect_error) {
    die("Connection failed: " . $conn->connect_error);
}

$quoted_names = array_map(function ($name) {
    return "'" . trim($name) . "'";
}, $st_names);
$name_list = implode(',', $quoted_names);

$sql = "SELECT `alias`, `main_id` from StarAliases where `alias` in ($name_list)";

$result = $conn->query($sql);
$data = $result->fetch_all();

// Return JSON
echo json_encode($data);

$conn->close();
?>

