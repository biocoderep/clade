from clade.validation.stage3_direction import direction_check


def test_positive_coefficient_is_correct_direction():
    result = direction_check(5.14)
    assert result["correct_direction"] is True
    assert "correct" in result["verdict"]


def test_negative_coefficient_is_wrong_direction():
    result = direction_check(-1.59)
    assert result["correct_direction"] is False
    assert "wrong" in result["verdict"]


def test_missing_coefficient_returns_none():
    result = direction_check(None)
    assert result["correct_direction"] is None
