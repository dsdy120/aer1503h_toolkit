import astropy
import astropy.coordinates
import astropy.time
import lamberthub
import numpy as np
import datetime
import orb_mech_constants

def planetary_position_velocity(iso_date:str, body:str):
    """
    Get the position of a celestial body at a given date.

    Parameters:
    date (datetime): The date for which to get the position.
    body (str): The name of the celestial body.

    Returns:
    numpy.ndarray: The position vector of the celestial body.
    """
    obs_time = astropy.time.Time(iso_date, format='isot', scale='utc')
    pos_vel = astropy.coordinates.get_body_barycentric_posvel(body,obs_time)
    pos = np.array([pos_vel[0].x.value, pos_vel[0].y.value, pos_vel[0].z.value])*orb_mech_constants.AU
    vel = np.array([pos_vel[1].x.value, pos_vel[1].y.value, pos_vel[1].z.value])*orb_mech_constants.AU/86400.0
    print(f"Position of {body:7s} at {iso_date} (JD: {obs_time.jd1} {obs_time.jd2: 5.1f}): {pos}, {vel}")

    return pos, vel

def lambert_interplanetary(start_date, end_date, body1, body2):
    """
    Calculate the Lambert trajectory between two celestial bodies.

    Parameters:
    start_date (datetime): The start date of the trajectory.
    end_date (datetime): The end date of the trajectory.
    body1 (str): The name of the first celestial body.
    body2 (str): The name of the second celestial body.

    Returns:
    tuple: A tuple containing the position and velocity vectors at the start and end dates.
    """
    time_of_flight = (datetime.datetime.fromisoformat(end_date) - datetime.datetime.fromisoformat(start_date)).total_seconds()
    departure_position,body1_velocity = planetary_position_velocity(start_date, body1)
    arrival_position,body2_velocity = planetary_position_velocity(end_date, body2)
    departure_velocity, arrival_velocity = lamberthub.gooding1990(orb_mech_constants.MU_SUN,departure_position, arrival_position, time_of_flight)
    return departure_position, departure_velocity, arrival_position, arrival_velocity, body1_velocity, body2_velocity


if __name__ == "__main__":
    # Example usage
    start_date = datetime.datetime(2005, 12, 1).isoformat()
    end_date = datetime.datetime(2006, 4, 1).isoformat()
    body1 = "Earth"
    body2 = "Venus"

    dep_pos, dep_vel, arr_pos, arr_vel, body1_vel, body2_vel = lambert_interplanetary(start_date, end_date, body1, body2)

    print("Departure Position:", dep_pos)
    print("Departure Velocity:", dep_vel)
    print("Arrival Position:", arr_pos)
    print("Arrival Velocity:", arr_vel)
    print("Body 1 Velocity:", body1_vel)
    print("Body 2 Velocity:", body2_vel)
    print("Departure V_inf:", np.linalg.norm(dep_vel - body1_vel))
    print("Arrival V_inf:", np.linalg.norm(arr_vel - body2_vel))