<?php 
header("Content-Type:application/json");
include("config.php"); ?>
<?php
$conn = new mysqli($servername, $username, $password, $dbname);
// Check connection
if ($conn->connect_error) {
    die("Connection failed: " . $conn->connect_error);
} 

$sql = 'SELECT 
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
    `st_radv`,
    `st_psfgrade_nfb1_high`,
    `st_psfgrade_nfb1_med`,
    `st_psfgrade_specb3_high`,
    `st_psfgrade_specb3_med`,
    `st_psfgrade_wfb4_high`,
    `st_psfgrade_wfb4_med`
    from Stars where sy_caltype = "RefStar"';
$result = $conn->query($sql);
$data = $result->fetch_all();

// Return JSON
echo json_encode($data);

$conn->close();
?>

