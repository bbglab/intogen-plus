process DriverSaturation {
	tag "Driver saturation"
	publishDir "${STEPS_FOLDER}/boostDM/saturation", mode: "copy"
	label "core"

    input:
        path (drivers) from DRIVERS_SUMMARY

    output:
		path("*.vep.gz") into DRIVERS_SATURATION

	script:
		"""

		drivers-saturation --drivers ${drivers}
		
		"""
}